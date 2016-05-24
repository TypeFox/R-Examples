#' @title Estimate adjacency matrix for equivalent FLC distributions based on states
#'
#' @description 
#' This function estimates the adjacency matrix \eqn{\mathbf{A}} of all pairwise
#' equivalent FLC distributions given the states \eqn{s_1, \ldots, s_K}. 
#' See Details below.
#' 
#' @section Details and user-defined distance function:
#' 
#' The \eqn{(i,j)}th element of the adjacency matrix is defined as
#' \deqn{
#'   \mathbf{A}_{ij} = distance(P(X \mid s_i), P(X \mid s_j)) = distance(f, g),
#' }
#' where \code{distance} is either 
#' \describe{
#'    \item{a metric}{in the function space of pdfs \eqn{f} and \eqn{g}, or}
#'    \item{a two sample test}{for \eqn{H_0: f=g}, e.g. a Kolmogorov-Smirnov 
#'                             test (\code{distance="KS"}).}
#' }
#' 
#' Again we use a functional programming approach and allow the user to specify 
#' any valid distance/similarity function \code{distance = function(f, g) 
#' return(...)}.
#' 
#' If \code{distance="KS"} the adjacency matrix contains p-values of a 
#' Kolmogorov-Smirnov test or the thresholded versions (if \code{alpha!=NULL})
#' - see Return for details.
#' 
#' Otherwise \code{distance} is an R function that takes as an input two vectors
#' \code{f} and \code{g} (e.g. the \code{\link{wKDE}} estimates for two 
#' states), and returns a non-negative, real number to estimate
#' their distance. Default is the \eqn{L_1} distance 
#' \code{distance = function(f, g) return(mean(abs(f-g)))}.
#' 
#' 
#' @param states vector of length \eqn{N} with entry \eqn{i} being the label 
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param FLCs \eqn{N \times n_f} matrix of FLCs (only necessary if 
#' \code{distance= "KS"})
#' @param pdfs.FLC \eqn{N \times K} matrix of all \eqn{K} state-conditional FLC 
#' densities evaluated at each FLC \eqn{\ell^{+}_i}, \eqn{i=1, \ldots, N} (only 
#' necessary if \code{distance = function(f, g) return(...)}).
#' @param alpha significance level for testing. Default: \code{alpha=NULL} 
#' (this will return a p-value matrix if \code{method == "KS"})
#' @param distance either a Kolmogorov-Smirnov test (\code{distance = "KS"}) 
#' or a function metric (e.g. \eqn{L_q} distance). For a distance function, 
#' \code{distance} requires as input a function of \eqn{f} and \eqn{g} that 
#' returns one value. 
#' 
#' Default: \code{distance = function(f, g) return(mean(abs(f-g)))} 
#' \eqn{\rightarrow} \eqn{L_1} distance.
#' @return
#' A \eqn{K \times K} adjacency matrix with a trimmed version of exp(-distance) or 
#' p-values.  If \code{alpha!=NULL} then it returns
#' the thresholded \eqn{0/1} matrix.  However, here \eqn{1} stands for equivalent, 
#' i.e. not rejecting. 
#' The matrix is obtained by checking for \code{pval>alpha} 
#' (rather than the usual \code{pval<alpha}).
#' @keywords manip nonparametric distribution multivariate
#' @export
#' @examples
#' WW = matrix(runif(10000), ncol = 10)
#' WW = normalize(WW)
#' temp_flcs = cbind(rnorm(nrow(WW)))
#' temp_pdfs.FLC = estimate_LC_pdfs(temp_flcs, WW)
#' AA_ks = estimate_state_adj_matrix(states = weight_matrix2states(WW), 
#'                                   FLCs = temp_flcs, distance = "KS")
#' AA_L1 = estimate_state_adj_matrix(pdfs.FLC = temp_pdfs.FLC)
#' 
#' par(mfrow = c(1,2), mar = c(1,1,2,1))
#' image2(AA_ks, zlim = c(0,1), legend = FALSE, main = "Kolmogorov-Smirnov")
#' image2(AA_L1, legend = FALSE, main = "L1 distance")
#' 

estimate_state_adj_matrix <- function(states = NULL, 
                                      FLCs = NULL, 
                                      pdfs.FLC = NULL, 
                                      alpha = NULL, 
                                      distance = function(f, g) return(mean(abs(f-g))) ) {
  if (is.character(distance)) {
    
    if (is.null(states) | is.null(FLCs)){
      stop("If you choose a KS test, then you must provide the state space
           assignment ('states') and the FLCs ('FLCs').")
    }
    
    num.states <- length(unique(states))
    adjacency.mat <- matrix(0, ncol = num.states, nrow = num.states)
    for (ii in seq_len(num.states)) {
      FLCs.in.ii <- FLCs[states == ii]
      for (jj in ii:num.states) {
        if (ii != jj) {
          adjacency.mat[ii, jj] <- suppressWarnings(ks.test(FLCs.in.ii, 
                                                            FLCs[states == jj])$p.value)
        }
      }
    }
  } else {
    num.states <- ncol(pdfs.FLC)
    adjacency.mat <- matrix(0, ncol = num.states, nrow = num.states)
    for (ii in seq_len(num.states)) {
      pdf.FLC.in.ii <- pdfs.FLC[, ii]
      for (jj in ii:num.states) {
        if (ii != jj) {
          adjacency.mat[ii, jj] <- distance(pdf.FLC.in.ii, pdfs.FLC[, jj])
        }
      }
    }
  }
  # fill up the whole matrix (since it's symmetric)
  adjacency.mat <- adjacency.mat + t(adjacency.mat)
  if (is.character(distance)) {
    diag(adjacency.mat) <- 1
  } else {
    adjacency.mat <- exp(-adjacency.mat)
    adjacency.mat <- adjacency.mat - min(adjacency.mat)
    adjacency.mat <- adjacency.mat / max(adjacency.mat)
  }
  
  if (!is.null(alpha) & is.character(distance)) {
    # merging matrix
    adjacency.mat <- (adjacency.mat > alpha)
  }
  
  invisible(adjacency.mat)
} 