eblup.mse.f.c3.asyvarcovarmat <- function(lme.obj, n.i,...){
    mat <- matrix(0,nrow=2,ncol=2)
    colnames(mat) <- c("v","e")
    rownames(mat) <- c("v","e")
    var.v <- as.numeric(VarCorr(lme.obj)[,1])[1]
    var.e <- as.numeric(VarCorr(lme.obj)[,1])[2]
    alpha.i <- var.e + n.i * var.v
        #vv = asy var of ranefs: 7.2.27
    mat[1,1] <- 0.5 * sum(n.i^2 * alpha.i^-2)
        #ee = asy var of resids: 7.2.28 - is sqrt(var.e)^-4 correct?
    mat[2,2] <- 0.5 * sum((n.i-1) * sqrt(var.e)^-4 + alpha.i^-2)
    mat[1,2] <- mat[2,1] <- 0.5 * sum(n.i * alpha.i^-2)
    #class(mat) <- "eblup.mse.f"
    return(mat)
}


