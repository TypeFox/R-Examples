#' Logistic Factor Analysis without C++ Dependency
#'
#' Estimate populatoin structure in genome-wide genotype matrices.
#'
#' It performs the logistic factor analysis, similar to \code{lfa} function in the lfa package.
#' This function works without C++ dependencies.
#' However, it would be much slower, does not include any other LFA-related functions, checks, and warnings.
#'
#' @param x a matrix with \code{m} loci (rows) and \code{n} observations (columns).
#' @param d a number of logistic factors.
#' @param ltrace a logical indicator as to whether to print the progress.
#'
#' @return \code{lfa.corpcor} returns a \code{n*d} matrix of \code{d} logistic factors.
#' The last column is always an intercept term.
#'
#' @importFrom corpcor fast.svd
#' @export lfa.corpcor
#' @seealso \link{jackstraw.LFA}
#'
#' @examples
#' set.seed(1234)
#' ## simulate genotype data from a logistic factor model
#' m=5000; n=100; pi0=.9
#' m0 = round(m*pi0)
#' m1 = m-round(m*pi0)
#' B = matrix(0, nrow=m, ncol=1)
#' B[1:m1,] = matrix(runif(m1*n, min=-.5, max=.5), nrow=m1, ncol=1)
#' L = matrix(rnorm(n), nrow=1, ncol=n)
#' BL = B %*% L
#' prob = exp(BL)/(1+exp(BL))
#'
#' dat = matrix(rbinom(m*n, 2, as.numeric(prob)), m, n)
#' out = lfa.corpcor(x=dat, d=2)
lfa.corpcor <- function(x, d, ltrace = FALSE){
    n.sv = d - 1
    m = nrow(x)
    n = ncol(x)

    norm_x = t(scale(t(x), scale=FALSE, center=TRUE))
    mysvd = fast.svd(norm_x)

    mean_x = apply(x,1,mean)
    sd_x = apply(x,1,sd)

    rm(norm_x)
    d = mysvd$d[1:n.sv]
    u = mysvd$u[,1:n.sv]
    v = mysvd$v[,1:n.sv]
    rm(mysvd)

    z = u %*% diag(d, n.sv, n.sv)  %*% t(v)

    z = (z*sd_x) + mean_x
    z = z/2
    rm(u); rm(d); rm(v)

    if(ltrace){
        print(summary(as.vector(z)))
    }

    #The .Call() is equivalent to the following lines of R code:
    #ind = as.logical(.Call("lfa_threshold", z, 1/(2*n)))
    zmin = apply(z, 1, min)
    zmax = apply(z, 1, max)
    ind  = (zmax<(1-2/n)) & (zmin>(2/n))
    z = z[ind,]
    z = log(z/(1-z))

    if(ltrace){
        print(dim(z))
    }

    norm_z = t(scale(t(z), scale=TRUE, center=TRUE))
    v= fast.svd(norm_z)$v[,1:n.sv]
    v = cbind(v,1)
    return(v)
}

#' Compute Deviance for Logistic Factors
#'
#' This function computes deviance between the full model and the null (intercept-only) model.
#' It uses built-in R functions, namely \code{glm}; slow but no C++ dependencies.
#' Make sure that \code{LFr1} and \code{LFr0} do not have intercept terms.
#'
#' @param dat a matrix with \code{m} rows and \code{n} columns.
#' @param LFr1 alternative logistic factors (an output from lfa or lfa.corpcor)
#' @param LFr0 null logistic factors (an output from lfa or lfa.corpcor)
#' @param p estimate p-values (by default, "FALSE")
#'
#' @return When {p=FALSE} (by default), \code{dev.R} returns a vector of \code{m} deviances.
#' @return When {p=TRUE}, a list consisting of
#' \item{dev}{the \code{m} deviances}
#' \item{p.value}{the \code{m} p-values based on a chisq distribution}
#'
#' @import stats
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
dev.R = function(dat, LFr1, LFr0=NULL, p=FALSE) {
    n = ncol(dat)
    # if(is.null(LFr0)) LFr0 = 1
    if(is.null(LFr0)) LFr0 = matrix(1, n, 1)

    obs1 = apply(dat, 1, function(x) glm(cbind(x, 2-x) ~ LFr1, family=binomial(link="logit"))$deviance)
    obs0 = apply(dat, 1, function(x) glm(cbind(x, 2-x) ~ LFr0, family=binomial(link="logit"))$deviance)

    dev = obs0 - obs1

    if(p==TRUE) {
        p.value = 1-pchisq(dev,ncol(LFr1)-ncol(LFr0))
        return(list(p.value=p.value, dev=dev))
    } else {
        return(dev)
    }
}
