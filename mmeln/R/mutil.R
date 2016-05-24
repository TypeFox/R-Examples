#### Fonctions utiles pour la librairie MMELN.
#### Fonction de densite pour une multinormale par exemple.
#### Prenant plusieurs type de localisation par contre, on ne
#### peut que definir une matrice de covariance.



#### Densite d'une multinormale avec parametre de localisation Mu et
#### matrice de covariance Sigma.


dmnorm=function(X,Mu,Sigma)
{
    if(class(X)=="data.frame")
    {
        X=matrix(X)
    }
    p=numeric(1)
    if( (is.vector(X) || (is.matrix(X) && dim(X)[1]==1) )
       && (is.vector(Mu) || (is.matrix(Mu) && dim(Mu)[1]==1)) )
    {
        p=length(Mu)
        X=as.vector(X)
        Mu=as.vector(Mu)
        if(sum(is.na(X))>0)
        {
            p=p-sum(is.na(X))
            X=na.exclude(X)
            X.na=na.action(X)
            return((2*pi)^(-p/2)/sqrt(det(Sigma[-X.na,-X.na]))*exp(-1/2*t(X-Mu[-X.na])%*%solve(Sigma[-X.na,-X.na])%*%(X-Mu[-X.na])))
        }
        return((2*pi)^(-p/2)/sqrt(det(Sigma))*exp(-1/2*t(X-Mu)%*%solve(Sigma)%*%(X-Mu)))
    }
    else if( (is.vector(X) || (is.matrix(X) && dim(X)[1]==1) )
            && (is.matrix(Mu) && dim(Mu)[1]>1) )
    {
        p=dim(Mu)[2]
        X=matrix(rep(as.vector(X),dim(Mu)[1]),ncol=p,byrow=TRUE)
    }
    else if( is.matrix(X) &
            (is.vector(Mu) || (is.matrix(Mu) && dim(Mu)[1]==1)))
    {
        p=length(Mu)
        N=dim(X)[1]
        Mu=matrix(rep(as.vector(Mu),N),N,byrow=TRUE)
    }
    else{p=dim(X)[2]}
    if(sum(is.na(X))>0)
    {
        result=numeric(dim(X)[1])
        for(i in 1:dim(X)[1])
        {
            if(sum(is.na(X[i,]))>0 & sum(is.na(X[i,]))!=p)
            {
                k=p-sum(is.na(X[i,]))
                Xi=na.exclude(X[i,])
                Xi.na=na.action(Xi)
                result[i]=(2*pi)^(-k/2)/sqrt(det(Sigma[-Xi.na,-Xi.na]))*exp(-1/2*t(Xi-Mu[1,-Xi.na])%*%solve(Sigma[-Xi.na,-Xi.na])%*%(Xi-Mu[1,-Xi.na]))
            }
            else
            {
                result[i]=(2*pi)^(-p/2)/sqrt(det(Sigma))*exp(-1/2*t(X[i,]-Mu[i,])%*%solve(Sigma)%*%(X[i,]-Mu[i,]))
            }
        }
        return(result)
    }
    return((2*pi)^(-p/2)/sqrt(det(Sigma))*exp(-1/2*apply(((X-Mu)%*%solve(Sigma)*(X-Mu)),1,sum)) )

}

plot.mmeln=function(x,...,main="",xlab="Temps",ylab="Y",col=1:x$G,leg=TRUE)
{
    X=x
    predict=matrix(0,ncol=X$G,nrow=X$p)
    for(i in 1:X$G)
    {
        predict[,i]=X$Xg[[i]]%*%X$param[[1]][[i]]
    }
    if(leg)
        oldpar=par(oma=c(2,0,0,0),las=1)
    matplot(predict,main=main,xlab=xlab,ylab=ylab,type="b",lty=1,col=col)
    if(leg)
    {
        txt=""
        if(X$G==1)
        {
            txt="Grp 1: 100%"
            mtext(txt, side = 1 ,outer = TRUE,  col = col,at =.11, adj=0 )
        }
        else
        {
            P=cbind(1,exp(X$Z%*%matrix(X$param$tau,ncol=X$G-1)))/apply(cbind(1,exp(X$Z%*%matrix(X$param$tau,ncol=X$G-1))),1,sum)
            txt=paste("Grp ",1:X$G,": ",round(apply(P,2,mean)*100,digits=1),"%",sep="")
            mtext(txt, side = 1 ,outer = TRUE,col = col,at=.15*(0:(X$G-1))+.11, adj=0)
        }

    }
    if(leg)
        par(oldpar)
}



#### Fonction pour estimer les parametres de la multinomiale.

multnm=function(Post,tau0,Z,g,iterlim=100,tol=1e-8)
{

    #### Fonction pour chaque itération.
    mn=function(P,tau0,Z,g)
    {
        tau=matrix(tau0,ncol=g-1)
        eta=Z%*%tau
        pi=cbind(1,exp(eta))/apply(cbind(1,exp(eta)),1,sum)
        var=pi*(1-pi)
        tau1=numeric()
        for(i in 1:(g-1))
        {
            tau1=c(tau1,as.vector(tau[,i])+solve(t(Z)%*%diag(as.vector(var[,i+1]))%*%Z)%*%t(Z)%*%(P[,i+1]-pi[,i+1]))
        }
        tau1
    }


    for(i in 1:iterlim)
    {
        tau1=mn(Post,tau0,Z,g)
        if(max(abs(tau1-tau0))<tol)
        {
            return(tau1)
        }
        tau0=tau1
    }
    stop(paste("Le nombre d'iterations maximales de",iterlim,"est depasse."))
}


##### Fonction pour le traitement des jeux de donnees incomplet.


##### Fonction pour faire une grande matrice inverse, dans le probleme de trouver la covariance
##### dans un jeu de donnee avec des donnees manquantes.


Xinv=function(Y,Sigma)
{
    p=dim(Y)[2]
    vYi=!is.na(Y[1,])
    bigSigma=solve(Sigma[vYi,vYi])
    for(i in 2:dim(Y)[1])
    {
        vYi=!is.na(Y[i,])
        bigSigma=rbind(cbind(bigSigma,matrix(0,dim(bigSigma)[1],sum(vYi)))
        ,cbind(matrix(0,sum(vYi),dim(bigSigma)[1]),solve(Sigma[vYi,vYi])))
    }
    bigSigma
}


##### Fonction pour calculer des covariances ponderees avec des
##### donnees manquantes selon la methode pairwise. (seulement pour calculer des bons starting values)

covNA.wt=function(Y,wt)
{
    p=dim(Y)[2]
    S=sqrt(apply(Y,2,
    function(x){xn=!is.na(x);
                (t(x[xn]-(sum(diag(wt[xn])%*%x[xn])/sum(wt[xn]))) %*% diag(wt[xn]) %*% (x[xn]-(sum(diag(wt[xn])%*%x[xn])/sum(wt[xn]))))/sum(wt[xn]) }))
    r=numeric()
    for(i in 1:p-1)
    {
        for(j in (i+1):p)
        {
            yn=apply(!is.na(cbind(Y[,i],Y[,j])),1,sum)==2
            r=c(r,(t(Y[yn,i]-sum(Y[yn,i]*wt[yn])/sum(wt[yn])) %*% diag(wt[yn]) %*% (Y[yn,j]-sum(Y[yn,j]*wt[yn])/sum(wt[yn])))/sum(wt[yn])/S[i]/S[j])
        }
    }
    R.l=matrix(0,p,p)
    R.l[lower.tri(R.l)]=r
    diag(S)%*%(R.l+diag(p)+t(R.l))%*%diag(S)
}


logit=function(xi)
{
    exp(xi)/(1+exp(xi))
}
