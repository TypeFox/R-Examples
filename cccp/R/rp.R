##
## Function for risk-parity optimization
## 
rp <- function(x0, P, mrc, optctrl = ctrl()){
    n <- nrow(P)
    m <- ncol(P)
    if(!identical(n, m)){
        stop("Matrix 'P' must be square.\n")
    }
    P <- 2 * P
    mrc <- as.matrix(mrc)
    x0 <- as.matrix(x0)
    if(!identical(nrow(mrc), n)){
        stop("Count of marginal risk contributions is not equal to the problem size.\n")
    }
    if(!identical(nrow(x0), n)){
        stop("Count of start values 'x0' is not equal to the problem size.\n")
    }
    if(!all.equal(sum(mrc), 1.0)){
        warning("Sum of marginal risk contributions does not equal one: Normalizing risk contributions.\n")
        mrc <- mrc / sum(mrc)
    }
    if(class(optctrl) != "Rcpp_CTRL"){
        stop("Provided argument for 'optctrl' is not a reference class object 'Rcpp_CTRL'.\n")
    }
    rpp(x0, P, mrc, optctrl)
}

