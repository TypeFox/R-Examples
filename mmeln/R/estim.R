##########################################
### Fonction pour estimer les parametres #
#########################################

### Prototype : estim(X,start,iterlim,tol)

estim=function(X,...)
{
    UseMethod("estim")
}

estim.mmeln=function(X,...,mu=NULL,tau=NULL,sigma=NULL,random.start=FALSE,iterlim=500,tol=1e-8)
{
    if(random.start | (is.null(mu) & is.null(tau) & is.null(sigma)))
    {
        if(X$pm>0)
            tau=rep(0,X$pm)
        else
            tau=NULL
        muT=apply(X$Y,2,function(x){mean(x,na.rm=TRUE)})
        ST=var(X$Y,use="pairwise.complete.obs")
        RT=diag(1/sqrt(diag(ST)))%*%ST%*%diag(1/sqrt(diag(ST)))
        sigma=list()
        if(X$cov=="UN"){ sigma[[1]]=c(sqrt(diag(ST)),RT[lower.tri(RT)]) }
        else if(X$cov=="CS" | X$cov=="AR1")
        {
            sigma[[1]]=c(sqrt(mean(diag(ST))),mean(ST[lower.tri(ST)])/mean(diag(ST)))
        }
        else if(X$cov=="UCS" | X$cov=="UAR1")
        {
            sigma[[1]]=c(sqrt(diag(ST)),mean(ST[lower.tri(ST)])/mean(diag(ST)))
        }
        else if(X$cov=="IND"){ sigma[[1]]=c(sqrt(mean(diag(ST)))) }
        else if(X$cov=="UIND"){ sigma[[1]]=c(sqrt(diag(ST)))}
        else{ stop("This type of covariance is not yet implemented") }
        mu=list()
        MU=matrix(rnorm(X$G*X$p),X$G,X$p)%*%chol(ST)+matrix(muT,X$G,X$p)
        for(g in 1:X$G)
        {
            mu[[g]]=solve(t(X$Xg[[g]])%*%solve(ST)%*%X$Xg[[g]])%*%t(X$Xg[[g]])%*%solve(ST)%*%MU[g,]
            if(!X$equalcov & g>1)
                sigma[[g]]=sigma[[1]]
        }
    }
    #### Verification des parametres de localisation.
    if(!is.list(mu) | is.null(mu) | length(mu)!=X$G)
    {
        stop("mu must be a list of length X$G")
    }
    for(g in 1:X$G)
    {
        if(dim(X$Xg[[g]])[2]!=length(mu[[g]]))
        {
            stop(paste("Wrong number of location (mu) parameters in group",g))
        }
    }
    #### Verification des parametres du melange.
    if(X$G==1 & !is.null(tau))
    {
        stop("tau must be set to NULL when only 1 group is estimated")
    }
    else if(X$G>1 & (is.null(tau) | !is.vector(tau) | length(tau)!=X$pm))
    {
        stop(paste("tau must be of length",X$pm))
    }
    #### Verification des parametres de covariance.
    if(!X$equalcov){
        if(is.null(sigma) | !is.list(sigma) | length(sigma)!=X$G)
        {
            stop("sigma must be a list of length X$G")
        }
        for(g in 1:X$G)
        {
            if(is.null(sigma[[g]]) | !is.vector(sigma[[g]]) | length(sigma[[g]])!=(X$pc/X$G))
            {
                stop(paste("Wrong number of covariance (sigma) parameters in group",g))
            }
        }

    }
    else
    {
        if(is.null(sigma) | !is.list(sigma) | length(sigma[[1]])!=X$pc)
        {
            stop(paste("Wrong number of covariance parameters, with equal covariance sigma must be of length 1"))
        }
    }
    if(X$cov=="IND" & X$equalcov==FALSE)
    {
        result=estimmmelnIND1(X,list(mu=mu,tau=tau,sigma=sigma),iterlim,tol);
        X$param=list(mu=result[[1]],tau=result[[2]],sigma=result[[3]])
        X$niter=result[[4]]
        X$H1=I.IND1(X)
        X$H2=IE.IND1(X)
        class(X)=c("mmelnSOL","mmeln")
        X
    }
    else if(X$cov=="CS" & X$equalcov==FALSE)
    {
        result=estimmmelnCS1(X,list(mu=mu,tau=tau,sigma=sigma),iterlim,tol);
        X$param=list(mu=result[[1]],tau=result[[2]],sigma=result[[3]])
        X$niter=result[[4]]
        X$H1=I.CS1(X)
        X$H2=IE.CS1(X)
        class(X)=c("mmelnSOL","mmeln")
        X
    }
    else
    {
        stop("This type of covariance is not yet implemented")
    }
}

#### Fonction d'affichage pour une solution d'un melange.
print.mmelnSOL=function(x,...,se.estim="MLR")
{
    X=x
    #### Fonction interne intermediaire.
    prnt.loc.mmeln=function(X)
    {

        loc=numeric()
        se.loc=numeric()
        Pst=post(X)
        rownames=NULL
        if(se.estim=="MLR")
            se.loc=sqrt(diag(solve(X$H1)%*%X$H2%*%solve(X$H1))[1:X$pl])
        else if(se.estim=="ML")
            se.loc=sqrt(diag(solve(X$H1))[1:X$pl])
        else if(se.estim=="ML.E")
            se.loc=sqrt(diag(solve(X$H2))[1:X$pl])
        else
            stop("Choice of estimator for SE are \"MLR\", \"ML\" or \"ML.E\"")
        for(g in 1:X$G)
        {
            loc=c(loc,X$param$mu[[g]])
            if(!is.null(dimnames(X$Xg[[g]])[[2]]))
            {
                rownames=c(rownames,paste("G",g,":",dimnames(X$Xg[[g]])[[2]],sep=""))
            }
            else
            {
                rownames=c(rownames,paste("Par",1:dim(X$Xg[[g]])[2],".G",g),sep="")
            }

        }
        z=loc/se.loc
        pz=round(pnorm(-abs(z))/2,digits=5)
        out.loc=data.frame(loc,se.loc,z,pz)
        names(out.loc)=c("Param","SE(Param)","Z","P(|Z|>z)")
        row.names(out.loc)=rownames
        out.loc
    }

    prnt.mel.mmeln=function(X)
    {
        if(X$G==1)
        {
            ### on ne fait rien car il n'y a qu'un seul groupe.
        }
        else
        {
            tau=matrix(X$param$tau,ncol=(X$G-1))
            mel=X$param$tau
            if(se.estim=="MLR")
                semel=sqrt(diag(solve(X$H1)%*%X$H2%*%solve(X$H1))[-(1:(X$pl+X$pc))])
            else if(se.estim=="ML")
                semel=sqrt(diag(solve(X$H1))[-(1:(X$pl+X$pc))])
            else if(se.estim=="ML.E")
                semel=sqrt(diag(solve(X$H2))[-(1:(X$pl+X$pc))])
            eta=X$Z%*%tau
            rownames=NULL
            pi=cbind(1,exp(eta))/apply(cbind(1,exp(eta)),1,sum)

            for(g in 2:X$G)
            {
                if(is.null(dimnames(X$Z)[[2]]))
                {
                    rownames=c(rownames,paste("G",g," vs G1:P",1:dim(X$Z)[2],sep=""))
                }
                else
                {
                    rownames=c(rownames,paste("G",g," vs G1:",dimnames(X$Z)[[2]],sep=""))
                }
            }
            z=mel/semel
            pz=round(pnorm(-abs(z))*2,digits=5)
            out.mel=data.frame(mel,semel,z,pz)
            names(out.mel)=c("Param","SE(Param)","Z","P(|Z|>z)")
            row.names(out.mel)=rownames
            return(out.mel)
        }
    }


    if(X$pm!=0)
    {
        theta=X$Z%*%matrix(X$param[[2]],ncol=X$G-1)
        P=cbind(1,exp(theta))/apply(cbind(1,exp(theta)),1,sum)
    }
    cat("Solution of mixture estimation :\n\nNombre d'iterations : ")
    cat(X$niter)
    if(se.estim=="MLR")
        cat("\n\nLocation parameters (Robust s.e. estimator) :\n")
    else if(se.estim=="ML")
        cat("\n\nLocation parameters (ML s.e. estimator) :\n")
    else
        cat("\n\nLocation parameters (Empirical ML s.e. estimator) :\n")
    print(prnt.loc.mmeln(X))
    if(X$pm!=0)
    {
        cat("\n\nMixture parameters :\n")
        print(prnt.mel.mmeln(X))
        cat("\n\nProportions of the mixture :\n")
        Prop=as.data.frame(matrix(round(apply(P,2,mean)*100,digits=1),ncol=X$G))
        dimnames(Prop)=list("Prop(%)",paste("Group",1:X$G))
        print(Prop)
    }
    else
        cat("\n\nProportion in mixture :\nProp(%) 100")
    cat("\n\nVariance-covariance type is ")
    cat(X$cov)
    cat(" :\n")
    if(X$equalcov)
    {
        cat("\n\nHomogeneous Covariance across groups:\n")
        S=cov.tsf(X$param[[3]][[1]],X$cov,X$p)
        dimnames(S)=list(paste("T",1:X$p,sep=""),paste("T",1:X$p,sep=""))
        print(S)
    }
    else
    {
        for(i in 1:X$G)
        {
            cat("\n\nCovariance matrix in group ")
            cat(i)
            cat("\n")
            S=cov.tsf(X$param[[3]][[i]],X$cov,X$p)
            dimnames(S)=list(paste("T",1:X$p,sep=""),paste("T",1:X$p,sep=""))
            print(S)
        }
    }
    cat("\n\nNumber of parameters : ")
    cat(X$pl+X$pm+X$pc)
    cat("\n\nlog(Likelihood) : ")
    cat(logLik(X))
    cat("\n")
}


cov.tsf=function(param,type,p)
{
    if(type=="UN")
    {
        R.lower=matrix(0,p,p)
        R.lower[lower.tri(R.lower)]=param[(p+1):length(param)]
        R=R.lower+diag(p)+t(R.lower)
        return(diag(param[1:p])%*%R%*%diag(param[1:p]))
    }
    else if(type=="CS")
    {
        R.lower=matrix(0,p,p)
        R.lower[lower.tri(R.lower)]=param[2]
        R=R.lower+diag(p)+t(R.lower)
        return(diag(rep(param[1],p))%*%R%*%diag(rep(param[1],p)))
    }
    else if(type=="UCS")
    {
        R.lower=matrix(0,p,p)
        R.lower[lower.tri(R.lower)]=param[p+1]
        R=R.lower+diag(p)+t(R.lower)
        return(diag(param[1:p])%*%R%*%diag(param[1:p]))
    }
    else if(type=="AR1")
    {

        R.lower=matrix(0,p,p)
        R.lower[lower.tri(R.lower)]=param[2]
        R=R.lower+diag(p)+t(R.lower)
        R=R^abs(matrix(1:p,p,p,byrow=TRUE)-matrix(1:p,p,p))
        return(diag(rep(param[1],p))%*%R%*%diag(rep(param[1],p)))
    }
    else if(type=="UAR1")
    {
        R.lower=matrix(0,p,p)
        R.lower[lower.tri(R.lower)]=param[p+1]
        R=R.lower+diag(p)+t(R.lower)
        R=R^abs(matrix(1:p,p,p,byrow=TRUE)-matrix(1:p,p,p))
        return(diag(param[1:p])%*%R%*%diag(param[1:p]))
    }
    else if(type=="IND")
    {
        return(diag(rep(param^2,p)))
    }
    else if(type=="UIND")
    {
        return(diag(param^2))
    }
    else
    {
        stop("Specified covariance is not yet implemented")
    }
}
