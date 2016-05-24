    #third component
eblup.mse.f.c3 <- function(lme.obj, asympt.var.covar#asymptotic variance-covariance matrix of variance components
                     , n.i,...){
    var.v <- as.numeric(VarCorr(lme.obj)[,1])[1]
    var.e <- as.numeric(VarCorr(lme.obj)[,1])[2]
        #approx covariances - is that right?
    ## V.bar.vv <- lme.obj$apVar[1,1]
    ## V.bar.ee <- lme.obj$apVar[2,2]
    ## V.bar.ve <- lme.obj$apVar[1,2]
        #inverted apprx var-cov matrix
    inv.var.covar <- solve(asympt.var.covar)
    V.bar.vv <- inv.var.covar[1,1]
    V.bar.ee <- inv.var.covar[2,2]
    V.bar.ve <- inv.var.covar[1,2]
        #7.2.23
    h <- var.e^2 * V.bar.vv + var.v^2 * V.bar.ee - 2 * var.e * var.v * V.bar.ve
        #7.2.22
    res <- n.i^-2 * (var.v + var.e/n.i)^-3 * h#is ^-3 correct? not ^-4? see 7.2.32
    #class(res) <- "eblup.mse.f"
    return(res)
}


