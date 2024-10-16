#Imports
library(tidyverse)
library(DESeq2)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x 1+m) tibble with a 'gene' column followed by
#' sample names as column names.
#'
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#'
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){
  result <- read_delim(filename)
  return(result)
}


#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns with sample names as column names.
#'
#' @return tibble: a (n x 1+m) tibble with a 'gene' column followed by m columns
#' of raw counts with genes that have zero variance across samples removed
#'
#' @note (g >= n)
#'
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
  result <- verse_counts %>% 
    # manipulate data to get column for rowwise variance
    select(!gene) %>% 
    mutate(row = row_number()) %>%
    pivot_longer(-row) %>%
    group_by(row) %>%
    summarize(var = var(value)) %>% # variance calculation
    bind_cols(verse_counts, .) %>% 
    # filter out zero variance rows
    filter(var!=0) %>% 
    # remove temp columns
    select(!c(row, var))
  return(result)
}


#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

timepoint_from_sample <- function(x) {
  result <- substring(x, 2,3) # take just the timepoint part
  return(result)
}

#' Grab sample replicate number from sample name
#'
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

sample_replicate <- function(x) {
  result <- str_split(x, '_', simplify = TRUE)[, 2] # split by _ and get second value for the replicate
  return(result)
}



#' Generate sample-level metadata from sample names.
#'
#' Will include columns named "sample", "timepoint", and "replicate" that store
#' sample names, sample time points, and sample replicate, respectively.
#'
#'
#' @param sample_names vector: character vector of length (_S_) consisting of sample
#' names from count data.
#'
#' @return tibble: a (_S_ x 3) tibble with column names "sample",
#' "timepoint", and "replicate". "sample" holds sample_names; "timepoint"
#' stores sample time points; and "replicate" stores sample replicate
#'
#' @note _S_ < m
#'
#' @example `meta <- meta_info_from_labels(colnames(count_data)[colnames(count_data)!='gene'])`

meta_info_from_labels <- function(sample_names) {
  result <- tibble(
    sample = sample_names,
    timepoint = timepoint_from_sample(sample_names),
    replicate = sample_replicate(sample_names)
  )
  return(result)
}





#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @return tibble of read totals from each sample. A tibble can be `(1 x _S_)` 
#' with sample names as columns names OR `(_S_ x 2)` with columns ("sample", "value")
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
  result <- count_data %>% select(-gene) %>% 
    summarise(across(everything(), sum)) %>% # sum the columns
    pivot_longer(cols = everything(), names_to = "sample", values_to = "value") # pivot to desired form
  return(result)
}



#' Normalize raw count data to counts per million WITH pseudocounts using the
#' following formula:
#'     count / (sample_library_size/10^6)
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m columns of cpm normalized read counts
#'
#' @examples
#' `normalize_by_cpm(count_data)`
normalize_by_cpm <- function(count_data) {
  # number of reads in each sample
  reads_cts <- get_library_size(count_data)
  df <- as.data.frame(count_data)
  
  # go through samples and library sizes
  for (i in 1:length(reads_cts$sample)) {
    sample_name <- reads_cts$sample[i]
    library_size <- reads_cts$value[i]
    
    # count / (library size) * 10^6
    df[[sample_name]] <- df[[sample_name]] / library_size * 10**6
  }
  
  return(as_tibble(df))
}


#' Normalize raw count matrix using DESeq2
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts

#' @param meta_data tibble: sample-level information tibble corresponding to the
#' count matrix columns
#'
#' @return tibble: DESeq2 normalized count matrix
#' @export
#'
#' @examples
#' `deseq_normalize(count_data, meta_data)`
deseq_normalize <- function(count_data, meta_data) {
  mat <- count_data %>% select(-gene) %>% as.matrix() # matrix for countData
  rownames(mat) <- count_data %>%  pull(gene)
  # making a summarized experiment object
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = meta_data,
                                design = ~ 1) 
  dds <- DESeq(dds)
  # normalize
  normalized_counts <- counts(dds, normalized=TRUE) %>% 
    as.data.frame %>% 
    rownames_to_column(var = 'gene') %>% 
    as_tibble()
  return(normalized_counts)
}



#' Perform and plot PCA using processed data.
#'
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs.
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(data, meta, title="") {
  # use plotPCA from DESeq2
  # set up dds
  mat <- data %>% select(-gene) %>% as.matrix() # matrix for countData
  rownames(mat) <- data %>%  pull(gene)
  # making a summarized experiment object
  dds <- DESeqDataSetFromMatrix(countData = round(mat), # DESeqDataSet needs countData to be non-negative integers
                                colData = meta,
                                design = ~ timepoint) 
  dds <- DESeq(dds)
  
  # DESeq2 vignette and use either the VST or RLog transformation on the raw counts
  rld <- rlog(dds, blind=FALSE)
  
  # PCA and plot
  pcas <- plotPCA(rld, intgroup = "timepoint", returnData = TRUE)
  percent_var <- round(100 * attr(pcas, "percentVar"))  # Get variance explained by each PC

  # Step 4: Create PCA plot using ggplot2
  g <- ggplot(pcas, aes(x = PC1, y = PC2, color = timepoint)) +
    geom_point() +
    labs(title = title,
         x = paste0("PC1: ", percent_var[1], "% variance"),
         y = paste0("PC2: ", percent_var[2], "% variance"))
  return(g)
}


#' Plot gene count distributions for each sample using boxplots.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#'
#' @example `plot_sample_distributions(data, scale_y_axis=TRUE, title='Raw Count Distributions')`

plot_sample_distributions <- function(data, scale_y_axis=FALSE, title="") {
  data <- data %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "counts")
  g <- ggplot(data, aes(x=fct_inorder(sample), y=counts, color=sample)) +
    geom_boxplot() + 
    labs(title = title,
         x ="sample", y = "counts")
  
  if (scale_y_axis==TRUE){
    g <- g + scale_y_continuous(trans='log10')
  }
  return(g)
}


#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#'
#' @example `plot_variance_vs_mean(data, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`

plot_variance_vs_mean <- function(data, scale_y_axis=FALSE, title="") {
  data <- data %>% rowwise() %>%
    mutate(mean_value = mean(c_across(-gene)), # mean for genes
           variance = var(c_across(-gene))) %>% # var for genes
    ungroup() %>% 
    arrange(mean_value) %>% # order by rank mean count
    mutate(rank = row_number()) # set ranks
  
  # scatter plot
  g <- ggplot(data, aes(x = rank, y = variance)) +
    geom_point(alpha = 0.5) + # transparency increased 
    geom_smooth(aes(y = variance, x = rank), span = .2, se = FALSE, color = "blue") + # smooth line
    labs(title = title,
         x ="Rank(Mean)", y = "Variance")
  
  if (scale_y_axis==TRUE){
    g <- g + scale_y_continuous(trans='log10')
  }
  return(g)
}

