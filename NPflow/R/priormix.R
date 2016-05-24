#'Construction of an Empirical based prior
#'
#'@param sDPMclust an object of class \code{summary.DPMMclust}
#'
#'@param nu0add an additionnal value integer added to hyperprior parameter nu
#'(increase to avoid non positive definite matrix sampling)
#'
#'@seealso summary.DPMMclust
#'
#'@export
#'
#'@examples
#' rm(list=ls())
#'
#' #Number of data
#' n <- 2000
#' set.seed(123)
#' #set.seed(4321)
#'
#'
#' d <- 2
#' ncl <- 4
#'
#' # Sample data
#'
#' sdev <- array(dim=c(d,d,ncl))
#'
#' xi <- matrix(nrow=d, ncol=ncl, c(-1.5, 1.5, 1.5, 1.5, 2, -2.5, -2.5, -3))
#' #xi <- matrix(nrow=d, ncol=ncl, c(-0.5, 0, 0.5, 0, 0.5, -1, -1, 1))
#' psi <- matrix(nrow=d, ncol=4, c(0.4, -0.6, 0.8, 0, 0.3, -0.7, -0.3, -0.8))
#' nu <- c(100,15,8,5)
#' p <- c(0.15, 0.05, 0.5, 0.3) # frequence des clusters
#' sdev[, ,1] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.3))
#' sdev[, ,2] <- matrix(nrow=d, ncol=d, c(0.1, 0, 0, 0.3))
#' sdev[, ,3] <- matrix(nrow=d, ncol=d, c(0.3, 0.15, 0.15, 0.3))
#' sdev[, ,4] <- .3*diag(2)
#'
#'
#' c <- rep(0,n)
#' w <- rep(1,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  w[k] <- rgamma(1, shape=nu[c[k]]/2, rate=nu[c[k]]/2)
#'  z[,k] <- xi[, c[k]] + psi[, c[k]]*rtruncnorm(n=1, a=0, b=Inf, mean=0, sd=1/sqrt(w[k])) +
#'                 (sdev[, , c[k]]/sqrt(w[k]))%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
#'  #cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'
#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["b_xi"]] <- rowMeans(z)
#' hyperG0[["b_psi"]] <- rep(0,d)
#' hyperG0[["kappa"]] <- 0.001
#' hyperG0[["D_xi"]] <- 100
#' hyperG0[["D_psi"]] <- 100
#' hyperG0[["nu"]] <- d+1
#' hyperG0[["lambda"]] <- diag(apply(z,MARGIN=1, FUN=var))/3
#'
#'  # hyperprior on the Scale parameter of DPM
#'  a <- 0.0001
#'  b <- 0.0001
#'
#'  nbclust_init <- 30
#'
#'\dontrun{
#'  MCMCsample_st <- DPMGibbsSkewT(z, hyperG0, a, b, N=2000, doPlot=FALSE,
#'                                 nbclust_init, diagVar=FALSE)
#'  s <- summary(MCMCsample_st, burnin = 1500, thin=5, posterior_approx=TRUE)
#'  pmix <- priormix(s)
#'  }
#'

priormix <- function(sDPMclust, nu0add=5){

    d <- nrow(sDPMclust$data)
    n <- ncol(sDPMclust$data)

    pmix <- list()
    oi <- sDPMclust$point_estim$opt_ind #optimal point estimate index
    pmix[["weights"]] <- as.vector(table(sDPMclust$point_estim$c_est)/n)
    USS <- sDPMclust$U_SS_list[[oi]]
    ncl <- length(USS)

    pmix[["alpha"]] <- ncl/log(n)

    pmix[["parameters"]] <- list()

    if(sDPMclust[["clust_distrib"]]=="skewT"){
        for(j in 1:ncl){
            pmix[["parameters"]][[j]] <- list()
            pmix[["parameters"]][[j]][["b_xi"]] <- USS[[j]][["xi"]]
            pmix[["parameters"]][[j]][["b_psi"]] <- USS[[j]][["psi"]]
            pmix[["parameters"]][[j]][["lambda"]] <- USS[[j]][["lambda"]]
            pmix[["parameters"]][[j]][["kappa"]] <- 0.001
            pmix[["parameters"]][[j]][["D_xi"]] <- 1
            pmix[["parameters"]][[j]][["D_psi"]] <- 1
            pmix[["parameters"]][[j]][["nu"]] <- d + 1 + nu0add
        }
    }
    return(pmix)
}