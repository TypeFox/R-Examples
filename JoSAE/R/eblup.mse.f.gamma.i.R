eblup.mse.f.gamma.i <-
function(lme.obj, n.i,...){
        #extract variance components: first for random, second for residual
    varcomp <- as.numeric(VarCorr(lme.obj)[,1])
    ret <- varcomp[1]/(varcomp[1]+varcomp[2]/n.i)#checked and ok
    #class(ret) <- "JoSAE.gamma.i.f"
    return(ret)
    }

