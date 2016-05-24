eblup.mse.f.c2.ai <- function(lme.obj, n.i, gamma.i, X.i,...)
{
        var.e <- as.numeric(VarCorr(lme.obj)[,1])[2]
        ## var.v <- as.numeric(VarCorr(lme.obj)[,1])[1]
        ## V.i <- var.e * diag(n.i) + var.v * rep(1,n.i) %*% t(rep(1,n.i))
        ## solve(V.i)#the same indeed
            #7.2.2
        V.i.inv <- 1/var.e * ( diag(n.i) - gamma.i/n.i * rep(1,n.i) %*% t(rep(1,n.i)))
            #7.2.7
        A.i <- t(X.i) %*% V.i.inv %*% X.i
        #class(A.i) <- "eblup.mse.f"
        return(A.i)
    }
