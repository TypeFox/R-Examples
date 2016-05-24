#################################################
##### Fonction pour calculer le log(Likelihood) #
#################################################


### Prototype : logLik.meln(X)

logLik.mmeln=function(object,...,param=NULL)
{
    X=object
    if(!is.null(param))
    {
        X$param=param
    }

    if(X$pm==0)
        P=rep(1,X$N)
    else
    {
        eta=X$Z%*%matrix(X$param[[2]],ncol=(X$G-1))
        P=cbind(1,exp(eta))/apply(cbind(1,exp(eta)),1,sum)
    }
    L=numeric()
    for( i in 1:X$G)
    {
        if(!X$equalcov)
            sig=cov.tsf(X$param[[3]][[i]],X$cov,X$p)
        else
            sig=cov.tsf(X$param[[3]][[1]],X$cov,X$p)
        L=cbind(L,dmnorm(X$Y,as.vector((X$Xg[[i]]%*%X$param[[1]][[i]])),sig))
        }
    lL=apply(P*L,1,function(x){log(sum(x))})
    return(sum(lL))
}

anova.mmeln=function(object,...,test=TRUE)
{
    X=object
    comp=list(...)
    sol=data.frame()
    if(test & length(comp)==1)
    {
        if(!inherits(comp[[1]],"mmeln"))
        {
            stop("Objects must be of class mmeln")
        }
        lL=c(logLik(X),logLik(comp[[1]]))
        df=c(X$pc+X$pm+X$pl,comp[[1]]$pc + comp[[1]]$pm + comp[[1]]$pl)
        N=c(sum(!is.na(X$Y)),sum(!is.na(comp[[1]]$Y)))
        N2=c(X$N,comp[[1]]$N)
        ent=c(entropy(X),entropy(comp[[1]]))
        sol=data.frame(df,lL,-2*lL,2*df-2*lL,-2*lL+df*log(N),-2*lL+df*log(N2),-2*lL+df*log(N)+2*ent
        ,c("",format(abs(diff(-2*lL)))),c("",format(1-pchisq(abs(diff(-2*lL)),abs(diff(df))))))
        names(sol)=c("df","logLik","-2*logLik","AIC","BIC","BIC2","ICL-BIC","L.Ratio","p-value")
        row.names(sol)=unlist(lapply(as.list(sys.call()[-1]),deparse))[1:dim(sol)[1]]
        class(sol)=c("anova.mmeln","data.frame")
        return(sol)
    }
    if(length(comp)==0)
    {
        lL=logLik(X)
        df=X$pc+X$pm+X$pl
        ent=entropy(X)
        N=sum(!is.na(X$Y))
        sol=data.frame(df,lL,-2*lL,2*df-2*lL,-2*lL+df*log(N),-2*lL+df*log(X$N),-2*lL+df*log(N)+2*ent)
        names(sol)=c("df","logLik","-2*logLik","AIC","BIC","BIC2","ICL-BIC")
        row.names(sol)=unlist(lapply(as.list(sys.call()[-1]),deparse))[1:dim(sol)[1]]
    }
    else
    {
        lL=df=N=N2=ent=vector("numeric",length(comp)+1)
        for(i in 1:(length(comp)+1))
        {
            if(i==1)
            {
                lL[i]=logLik(X)
                df[i]=X$pc+X$pm+X$pl
                N[i]=sum(!is.na(X$Y))
                N2[i]=X$N
                ent[i]=entropy(X)
            }
            else
            {
                if(!inherits(comp[[i-1]],"mmeln"))
                    stop("Objects must be of class mmeln")
                lL[i]=logLik(comp[[i-1]])
                df[i]=comp[[i-1]]$pc+comp[[i-1]]$pm+comp[[i-1]]$pl
                N[i]=sum(!is.na(comp[[i-1]]$Y))
                N2[i]=comp[[i-1]]$N
                ent[i]=entropy(comp[[i-1]])
            }
        }
        sol=data.frame(df,lL,-2*lL,2*df-2*lL,-2*lL+df*log(N),-2*lL+df*log(N2),-2*lL+df*log(N)+2*ent)
        names(sol)=c("df","logLik","-2*logLik","AIC","BIC","BIC2","ICL-BIC")
        row.names(sol)=unlist(lapply(as.list(sys.call()[-1]),deparse))[1:dim(sol)[1]]

    }
    sol
}
