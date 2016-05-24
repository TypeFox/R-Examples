    #the star-third componetnt
eblup.mse.f.c3.star <- function(lme.obj, asympt.var.covar#=lme.obj$apVar
                     , n.i, mean.resid.i,...){
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
        #7.2.23 this is how Rao suggests it - correct?
    h <- var.e^2 * V.bar.vv + var.v^2 * V.bar.ee - 2 * var.e * var.v * V.bar.ve
        #is it not correct to multiply the vars with the asympt vars?
    #h <- var.e^2 * V.bar.ee + var.v^2 * V.bar.vv - 2 * var.e * var.v * V.bar.ve
        #7.2.32
    res <- n.i^-2 * (var.v + var.e/n.i)^-4 * h * mean.resid.i^2#what is correct? see 7.2.22
    #class(res) <- "eblup.mse.f"
    return(res)
}

