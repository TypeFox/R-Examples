#' Determines species turnover
#'
#' 'Turnover' is a function that assesses species turnover from the range
#' perspective (traditional method).
#'
#' If the 'community' perspective is desired, simply transpose the matrix
#' before analysis using the transpose function ('t()'), but make sure you
#' understand the implications of this action, as the interpretation of the
#' output changes dramatically.
#'
#' 'method' is an argument handed to functions in the 'vegan' package. Leibold
#' & Mikkelson advocated the use of equiprobable rows and columns (provided
#' that rows and columns had at least one entry). This method is called 'r00'.
#' Methods maintaining row (site) frequencies include 'r0','r1' & 'r2'. The
#' default method argument is 'r1', which maintains the species richness of a
#' site (row totals) and fills species ranges (columns) based on their marginal
#' probabilities. Arguably the most conservative null algorithm is the fixed
#' row - fixed column total null, which is implemented as 'fixedfixed'. See the
#' help file for 'commsimulator' or Wright et al. 1998 for more information.
#'
#' If 'order' is FALSE, the interaction matrix is not ordinated, allowing the
#' user to order the matrix based on site characteristics or other biologically
#' relevant characteristics.
#'
#' This function can either be used as a standalone, or can be used through the
#' 'metacommunity()' function, which determines all 3 elements of metacommunity
#' structure (coherence, boundary clumping, & turnover) (Leibold & Mikkelson
#' 2002). The turnover metric used here is equivalent to the number of checkerboard
#' units community with species ranges (range perspective) filled in
#'
#' @param comm community data in the form of a presence absence matrix
#' @param method null model randomization method used by 'nullmaker'. See
#' details.
#' @param sims number of simulated null matrices to use in analysis
#' @param scores axis scores to ordinate matrix. 1: primary axis scores
#' (default) 2: secondary axis scores
#' @param order logical argument indicating whether to ordinate the interaction
#' matrix or not. See details.
#' @param allowEmpty logical argument indicating whether to allow null
#' matrices to have empty rows or columns
#' @param binary logical argument indicating whether to ordinate the community
#' matrix based on abundance or binary (default) data.
#' @param verbose Logical. Prints a graphical progress bar that tracks the
#' creation of null matrices. Useful for conservative null models on large
#' and/or sparse data.
#' @return A data.frame containing the test statistic (turnover), z-value (z), p-value (pval), mean (simulatedMean) and variance (simulatedVariance) of simulations, and randomization method (method)
#'
#' @author Tad Dallas
#' @export
#' @references Leibold, M.A. and G.M. Mikkelson. 2002. Coherence, species
#' turnover, and boundary clumping: elements of meta-community structure. Oikos
#' 97: 237 - 250.
#'
#' Wright, D.H., Patterson, B.D., Mikkelson, G.M., Cutler, A. & Atmar, W.
#' (1998). A comparative analysis of nested subset patterns of species
#' composition. Oecologia 113, 1-20.
#' @examples
#'
#' #define an interaction matrix
#' data(TestMatrices)
#' intmat <- TestMatrices[[3]]
#'
#' #determine species turnover
#' turnover.intmat <- Turnover(intmat, method='r1', sims=100, scores=1, binary=TRUE)
#'

Turnover <-function(comm ,method="r1" ,sims=1000 ,scores=1, order=TRUE, allowEmpty=FALSE, binary=TRUE, verbose=FALSE){

 if(order){comm = OrderMatrix(comm, scores = scores, binary = binary)}
 for(i in 1:ncol(comm)) {
   comm[min(which(comm[,i] == 1)):max(which(comm[,i] == 1)), i] <- 1
 }

 turnover <- function(web){
  D <- designdist(web, method = "(A-J)*(B-J)", terms = "minimum")
  return(sum(D))
 }

 statistic <- turnover(comm)
 nulls <- NullMaker(comm = comm, sims = sims, method = method,
        allowEmpty = allowEmpty, verbose = verbose, ordinate=order)
 simstat <- as.numeric(lapply(nulls,turnover))
 varstat <- sd(simstat)
 z <- (mean(simstat)-statistic)/(varstat)
 pval <- 2*pnorm(-abs(z))
 return(data.frame(turnover=statistic, z=z,pval=pval, simulatedMean=mean(simstat), simulatedVariance=varstat, method=method))
}
