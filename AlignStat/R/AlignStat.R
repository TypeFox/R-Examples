#' AlignStat: A tool for the statistical comparison of alternative multiple sequence alignments
#' 
#' This package contains functions that compare two alternative multiple sequence 
#' alignments (MSAs) to determine whether they align homologous residues in the same
#' columns as one another. It then classifies similarities and differences into 
#' conserved gaps, conserved sequence, merges, splits or shifts of one MSA relative 
#' to the other. Summarising these categories for each MSA column yields information
#' on which sequence regions are agreed upon my both MSAs, and which differ. Several
#' plotting functions enable easily visualisation of the comparison data for analysis.
#' 
#' @section Computing statistics:
#' Use \code{\link{compare_alignments}} to calculate statistics comparing two MSAs
#'
#' @section Plotting functions:
#' Use \code{\link{plot_similarity_heatmap}} to view a heatmap comparing the similarity
#' Use \code{\link{plot_dissimilarity_matrix}} to view a matrix of dissimilarity types
#' Use \code{\link{plot_similarity_summary}} to view a summary of column average similarity 
#' Use \code{\link{plot_dissimilarity_summary}} to view a summary of column dissimilarity types
#' 
#' 
#' @docType package
#' @name AlignStat
#' @useDynLib AlignStat
#' @importFrom Rcpp sourceCpp
NULL
