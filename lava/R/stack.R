##' Stack estimating equations
##'
##' Stack estimating equations
##' @param x Model 1
##' @param model2 Model 2
##' @param D1u Derivative of score of model 2 w.r.t. parameter vector of model 1
##' @param inv.D2u Inverse of deri
##' @param weight weight (vector or function)
##' @param dweight derivative of weight wrt parameters of model 1
##' @param U Optional score function (model 2) as function of all parameters
##' @param k Debug argument
##' @param keep1 If TRUE only parameters of model 2 i s returned
##' @param ... Additional arguments to lower level functions
##' @aliases stack.estimate
##' @export
stack.estimate <- function(x,model2,D1u,inv.D2u,weight,dweight,U,k=1,keep1=FALSE,...) {
    iid1 <- iid(x)
    iid2 <- iid(model2)
    if (missing(inv.D2u)) {
        inv.D2u <- -attributes(iid2)$bread
    }
    if (is.null(inv.D2u)) stop("Need derivative of second stage score")
    if (!missing(U)) {
        D1u <- numDeriv::jacobian(U,coef(x))
    }
    if (!missing(weight) && is.function(weight)) {
        dweight <- numDeriv::jacobian(weight,coef(x))
        weight <- weight(coef(x))
    }
    if (!missing(dweight)) {
        D2u <- Inverse(inv.D2u)        
        u2 <- iid2%*%D2u ## Score of stage two equation derived from estimated influence function
        ## Derivative of score wrt first set of parameters (weight-parameters)
        D1u <- crossprod(apply(u2,2,function(x) -x/weight),dweight)
    }
    ii <- iid(merge(x,model2))
    iid1. <- ii[,seq_along(coef(x)),drop=FALSE]
    iid2. <- ii[,length(coef(x))+seq_along(coef(model2)),drop=FALSE]
    iid3 <- t(inv.D2u%*%(D1u%*%t(iid1.)))
    if (!keep1) return(estimate(coef=coef(model2),iid=cbind(iid2.+k*iid3)))
    estimate(coef=c(coef(x),coef(model2)),iid=cbind(iid1.,iid2. + k*iid3))
}

##' @export
stack.glm <- function(x,model2,...) {
    stack(estimate(x),estimate(model2),...)
}


