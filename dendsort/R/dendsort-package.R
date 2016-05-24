#' Modular Leaf Ordering Methods for Dendrogram Nodes
#'
#' This package includes functions to optimize ordering of nodes in a dendrogram, without affecting the meaning of the dendrogram. 
#' A dendrogram can be sorted based on the average distance of subtrees, or based on the smallest distance value. These sorting 
#' methods improve readability and interpretability of tree structure, especially for tasks such as comparison of different 
#' distance measures or linkage types and identification of tight clusters and outliers. As a result, it also introduces more 
#' meaningful reordering for a coupled heatmap visualization.
#'
#'
#' @name dendsort-package
#' @aliases dendsort-package
#' @docType package
#' @title Modular Leaf Ordering Methods for Dendrogram Nodes
#' @importFrom stats as.dendrogram as.hclust is.leaf
#' @author Ryo Sakai \email{ryo.sakai@esat.kuleuven.be}
#' @keywords package
NULL


#' Sample data matrix from the integrated pathway analysis of gastric cancer from the Cancer Genome Atlas (TCGA) study
#'
#' a multivariate table obtained from the integrated pathway analysis of gastric cancer from the Cancer Genome Atlas (TCGA) study.
#' In this data set, each column represents a pathway consisting of a set of genes and each row represents a cohort of samples based
#' on specific clinical or genetic features. For each pair of a pathway and a feature, a continuous value of between 1 and -1 is
#' assigned to score positive or negative association, respectively.
#'
#' We would like to thank Sheila Reynolds and Vesteinn Thorsson from the Institute for Systems Biology for sharing this sample data set.
#'
#' @docType data
#' @keywords datasets
#' @name sample_tcga
#' @usage data(sample_tcga)
#' @format A data frame with 215 rows and 117 variables
NULL

