"periodhist" <- function(X,dfreq=FALSE,vt,drop=TRUE)
{

    ############################################
    # Validation des arguments fournis en entrée
    valid.one(dfreq,"logical")
    valid.vt(vt)
    Xvalid <- valid.X(X=X, dfreq=dfreq, vt=vt)
        X <- Xvalid$X
    valid.one(drop,"logical")
    ############################################

        I <- length(vt)
        Xinter <- if(dfreq) matrix(0,dim(X)[1],I+1) else matrix(0,dim(X)[1],I)
        for (i in (1:I))
        {
                if (i==1) { Xs <- matrix(X[,c(1:vt[i])],nrow=dim(X)[1]) } else
                { Xs <- matrix(X[,c((sum(na.rm=TRUE,vt[1:(i-1)])+1):sum(na.rm=TRUE,vt[1:i]))],nrow=dim(X)[1]) }
                Xinter[,i] <- apply(Xs,1,max)
        }
        if (dfreq) Xinter[,I+1] <- X[,sum(na.rm=TRUE,vt)+1]
        Y <- histfreq.t(Xinter,dfreq=dfreq)
        Xinterfreq <- cbind(histpos.t(I),Y)
        if (drop) Xinterfreq <- Xinterfreq[Y!=0,]
        return(Xinterfreq)
}
