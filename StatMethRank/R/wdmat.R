#' Compute the weighted distances between two data sets of rankings
#'
#' This function computes the weighted distances between two data sets 
#' of complete rankings. The results are put in the matrix form. The 
#' data set could be aggregated or not.
#'
#' @param dset1 one data set with each row being a ranking
#' @param dset2 the other data set (default as dset1)
#' @param dset1.agg whether the data set is in the aggregated form
#'                  (default as FALSE)
#' @param dset2.agg whether the data set is in the aggregated form
#'                  (default as FALSE)
#' @param dtype type of the weighted distance measure 
#' Kendall or K(default) : "Weighted Kendall's tau", SqrtSpearman 
#' or SS : "Square root of weighted Spearman", Spearman or S : 
#' "Weighted Spearman"
#' @param weight weight vector (default as all components being 1)
#' @param modal the modal ranking 
#'              (default as c(1:k), k being the number of ranked items)
#' @return a list whose first object is a vector about the aggregation 
#'         status of the two data sets c(dset1.agg, dset2.agg), and 
#'         second object is a matrix of distances.
#' @export
#' @author Yumin Zhang <zymneo@@gmail.com>
#' @examples
#' data(Croon)
#' wdmat(Croon, dset1.agg=TRUE, dset2.agg=TRUE)
wdmat <- function(dset1, dset2=dset1, dset1.agg=FALSE, dset2.agg=FALSE,
                  dtype="Kendall", weight=NULL, modal=NULL) 
{
	# Compute the weighted distances between two data sets of rankings
	#
	# This function computes the weighted distances between two data sets 
	# of complete rankings. The results are put in the matrix form. The 
	# data set could be aggregated or not.
	# Author 	: Yumin Zhang
	# Email 	: zymneo@gmail.com
	# Created 	: Oct 16, 2013
    if (is.null(dim(dset1))) {
        dset1 <- t(as.matrix(dset1))
    }
    if (is.null(dim(dset2))) {
        dset2 <- t(as.matrix(dset2))
    }
    # convert dype of dset to matrix (if, say, dset is of type array)
    if (!is.matrix(dset1)) {
        dset <- as.matrix(dset1)
    }
    if (!is.matrix(dset2)) {
        dset <- as.matrix(dset2)
    }
    # find nitem (integer, number of ranked items) 
    # assuming this number for dset1 & dset2 are the same
    if (dset1.agg) {
        nitem <- ncol(dset1) - 1
        dset1.original <- dset1
        dset1 <- dset1[, 1:nitem]  # (matrix) set of unique observed rankings
        freq1 <- dset1[, ncol(dset1)]  # (vector) frequencies of 
        # unique observed rankings
    } else {
        nitem <- ncol(dset1)
    }
    if (dset2.agg) {
        dset2.original <- dset2
        dset2 <- dset2[, 1:nitem]  # (matrix) set of unique observed rankings
        freq2 <- dset2[, ncol(dset2)]  # (vector) frequencies of 
        # unique observed rankings
    }
    # initialize weight and modal if they are not specified
    if (is.null(weight)) {
        weight <- rep(1, nitem)
    }
    if (is.null(modal)) {
        modal <- 1:nitem
    }
    dmat.row <- nrow(dset1)
    dmat.col <- nrow(dset2)
    dmat.result <- list(agg=c(dset1.agg, dset2.agg), 
                        mat=matrix(0, nrow=dmat.row, ncol=dmat.col))
    #
    # Define function for computing pairwise distance 
    #
    # Square root of weighted Spearman
    if (dtype == "SqrtSpearman" | dtype == "SS") {
        pair.dist <- function(x, y, weight, modal) {
            d <- 0
            for (j in 1:nitem) {
                d <- d + (x[j] - y[j])^2 * weight[modal[j]]
            }
            return(d^0.5)
        }
    } 
    # Weighted Spearman
    if (dtype == "Spearman" | dtype == "S") {
        pair.dist <- function(x, y, weight, modal) {
            d <- 0
            for (j in 1:nitem) {
                d <- d + (x[j] - y[j])^2 * weight[modal[j]]
            }
            return(d)
        }
    }
    # Weighted Spearman's footrule
    if (dtype == "Footrule" | dtype == "F") {
        pair.dist <- function(x, y, weight, modal) {
            d <- 0
            for (j in 1:nitem) {
                d <- d + abs(x[j] - y[j]) * as.numeric(weight[modal[j]]) 
            }
            return(d)
        }
    }
    # Weighted Kendall's tau
    if (dtype == "Kendall" | dtype == "K") {
        pair.dist <- function(x, y, weight, modal) {
            d <- 0
            for (j in 1:(nitem - 1)) {
                for (k in (j + 1):nitem) {
                    if ((x[j] - x[k]) * (y[j] - y[k]) < 0) {
                        d <- d + weight[modal[k]] * weight[modal[j]]
                    }
                }
            }
            return(d)
        }
    }
    for (j in 1:dmat.col) {
        dmat.result$mat[, j] <- dmat.result$mat[, j] + 
            t(apply(dset1, 1, pair.dist, y=dset2[j, ], 
                    weight=weight, modal=modal))
    }
    return(dmat.result)
}
