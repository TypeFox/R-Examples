#'ELoss of a partition point estimate compared to a gold standard
#'
#'Evaluate the loss of a point estimate of the partition compared to a gold standard according to a given loss function
#'
#'@param c vector of length \code{n} containing the estimated partition
#'of the  \code{n} observations.
#'
#'@param gs vector of length \code{n} containing the  gold standard
#'partition of the  \code{n} observations.
#'
#'@param lossFn character string specifying the loss function to be used.
#'Either "F-measure" or "Binder" (see Details). Default is "F-measure".
#'
#'@param a only relevant if \code{lossFn} is "Binder". Penalty for wrong
#'coclustering in \code{c} compared to code{gs}. Defaults is 1.
#'
#'@param b only relevant if \code{lossFn} is "Binder". Penalty for missed
#'coclustering in \code{c} compared to code{gs}. Defaults is 1.
#'
#'@return the cost of the point estimate \code{c} in regard of the
#'gold standard \code{gs} for a given loss function.
#'
#'@details The cost of a point estimate partition is calculated using either a pairwise
#' coincidence loss function (Binder), or 1-Fmeasure (F-measure).
#'
#'
#'@author Boris Hejblum
#'
#'@export
#'
#'@references J.W. Lau & P.J. Green. Bayesian Model-Based Clustering
#'Procedures, Journal of Computational and Graphical Statistics,
#'16(3): 526-558, 2007.
#'
#' D. B. Dahl. Model-Based Clustering for Expression Data via a
#' Dirichlet Process Mixture Model, in Bayesian Inference for
#' Gene Expression and Proteomics, K.-A. Do, P. Muller, M. Vannucci
#' (Eds.), Cambridge University Press, 2006.
#'
#'@seealso \code{\link{similarityMat}}, \code{\link{cluster_est_binder}}
#'

evalClustLoss <- function(c, gs, lossFn="F-measure", a=1, b=1){
    n <- length(c)
    loss <- NA

    if(length(gs)!=n){
        stop("'c' and 'gs' arguments have not the same length")
    }

    if(lossFn == "F-measure"){
        loss <- 1 - FmeasureC(pred=c, ref=gs)
    }else if(lossFn == "Binder"){
        c_coclust <- sapply(c, FUN=function(x){x==c})
        gs_coclust <- sapply(gs, FUN=function(x){x==gs})

        dif <- c_coclust-gs_coclust
        dif[which(dif==1)] <- b
        dif[which(dif==-1)] <- a

        loss <- sum(dif)

    }else{
        stop("Specified loss function not available.\n
             Specify either 'F-measure' or 'Binder' for the lossFn argument.")
    }

    return(loss)
}