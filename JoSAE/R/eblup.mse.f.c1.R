eblup.mse.f.c1 <-
function(lme.obj, n.i, gamma.i, ...){
    var.e <- as.numeric(VarCorr(lme.obj)[,1])[2]
    ret <- gamma.i * (var.e / n.i)
    #class(ret) <- "eblup.mse.f"
    return(ret)
}

