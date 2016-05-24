#' Functions for drawing gapped cluster heatmap with ggplot2
#'
#' This is a set of tools for drawing gapmaps using \code{\link[ggplot2]{ggplot}}
#' 
#' \code{\link{gap_data}} extracts data from a dendrogram object. Make sure to convert \code{hclust} object to \code{dendrogram} object by calling \code{as.dendrogram()}.
#' This method generates an object class \code{gapdata}, consisting of a list of \code{data.frames}.
#' The general workflow is as following:
#' \enumerate{
#' \item{Hierarchical clustering \code{hclust()}}
#' \item{Convert the \code{hclust} output class into \code{dendrogram} by calling \code{as.dendrogram()}}
#' \item{Generate a gapped cluster heatmap by specifying a \code{matrix} and \code{dendrogram} objects for rows and columns in \code{gapmap()} function}
#' }
#' 
#' 
#' @name gapmap-package
#' @aliases gapmap-package
#' @docType package
#' @title Draws gapped heatmap (gapmap) and gapped dendrograms using ggplot2 in [R].
#' @importFrom ggplot2 ggplot geom_segment theme element_blank aes_string coord_flip scale_x_continuous scale_y_reverse scale_y_continuous element_text element_line labs geom_tile scale_fill_gradientn geom_text
#' @importFrom grid unit arrow 
#' @importFrom stats is.leaf
#' @importFrom reshape2 melt
#' @author Ryo Sakai \email{ryo.sakai@esat.kuleuven.be}
#' @keywords package
NULL

#' Sample data matrix from the integrated pathway analysis of gastric cancer from the Cancer Genome Atlas (TCGA) study.
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