#################################################
### Fonctions reliee a l'estimation de modele ###
### avec une covariance de type IND.         ###
###                                         ###
### Analyste: Charles-Edouard Giguere      ###
#############################################


##### Algorithme EM applique a des melanges multinormaux avec covariance IND et covariance inegale.


estimmmelnIND1=function(X,param,iterlim,tol)
{
##### on calcule la fonction de log-vraissemblance initiale.
    logL0=logLik(X,param=param)
    logL1=mu1=sigma1=tau1=numeric()
    mu1=param[[1]]
    tau1=param[[2]]
    sigma1=param[[3]]
    for(iterationEM in 1:iterlim)
    {

##### Etape E : on calcule l'esperance d'etre dans l'un ou l'autre des groupes ce qui donne
#####           la probabilite a posteriori qui est utilise pour construire la fonction objective
#####           Q* a l'etape M.
        Post=post(X,tau=tau1,mu=mu1,sigma=sigma1)

##### Etape M : Calcul des proportions

        if(X$G>1)
        {
##### La fonction multnm se trouve dans le fichier mutil.R
            tau1=multnm(Post,tau1,X$Z,X$G,iterlim=iterlim,tol=tol)
        }
        else
        {
            tau1=NULL
        }
#### Etape M : Calcul des parametres de localisation et de dispersion.
        loc.disp=estimloc.disp.IND1(X,Post,iterlim,tol,iterationEM)
        mu1=loc.disp[[1]]
        sigma1=loc.disp[[2]]
        logL1=logLik(X,param=list(mu1,tau1,sigma1))
        if(X$G==1)
        {
#### on retourne param
            return(list(mu=mu1,tau=tau1,sigma=sigma1,iterationEM))
        }
        else if(abs((logL1-logL0)/logL0) < tol )
        {
#### on retourne param
            return(list(mu=mu1,tau=tau1,sigma=sigma1,iterationEM))
        }

        else if(is.infinite(logL0))
        {
            logL0=logL1;
        }
        else
        {
            logL0=logL1;
        }
    }
    warning("Algorithm reached iteration limit without convergence; Results may not be reliable")
    return(list(mu=mu1,tau=tau1,sigma=sigma1,iterlim))
}


#### Estimation du beta et de sigma.

estimloc.disp.IND1=function(X,Post,iterlim,tol,iterEM)
{
    sigma=beta=list()
    for(g in 1:X$G)
    {
        S1=S2=0
        for(i in 1:X$N)
        {
            Xij=X$Xg[[g]][X$Yv[[i]],]
            S1=S1+Post[i,g]*t(Xij)%*%Xij
            S2=S2+c(Post[i,g]*t(Xij)%*%X$Yl[[i]])
        }
        beta[[g]]=solve(S1,S2)
        S3=0
        for(i in 1:X$N)
        {
            Xij=as.matrix(X$Xg[[g]][X$Yv[[i]],])
            rij=(X$Yl[[i]]-Xij%*%beta[[g]])
            S3=S3+Post[i,g]*c(t(rij)%*%rij)
        }
        sigma[[g]]=sqrt(S3/c(t(X$pi)%*%Post[,g]))
    }
    return(list(beta,sigma))
}



#### methode qui calcule la matrice hessienne des parametres du melange.
#### idealement ne rouler qu'une fois car cette methode peut s'averer un
#### peu lente.

I.IND1=function(X)
{
    Hessian=matrix(0,X$pl+X$pm+X$pc,X$pl+X$pm+X$pc)
    Pst=post(X)
    if(X$G>1)
    {
        eta=X$Z%*%cbind(0,matrix(X$param$tau,ncol=X$G-1))
        P=exp(eta)/apply(exp(eta),1,sum)
    }
    beta=xi=lambda=beta.xi=beta.lambda=xi.lambda=list()
    for(i in 1:X$N)
        for(j in 1:X$G)
            for(k in 1:X$G)
                beta[[paste(j,k)]] = xi[[paste(j,k)]] = lambda[[paste(j,k)]] =
                    beta.xi[[paste(j,k)]] =  beta.lambda[[paste(j,k)]] = xi.lambda[[paste(j,k)]]=0
    for(i in 1:X$N)
    {
        pi=X$pi[i]
        Zi=X$Z[i,]
        for(j in 1:X$G)
            for(k in 1:X$G)
            {
                Xij=as.matrix(X$Xg[[j]][X$Yv[[i]],])
                rij=c(X$Yl[[i]]- c(Xij%*%X$param$mu[[j]]))
                sigj=X$param$sigma[[j]][1]
                if(j==k)
                {
                    beta[[paste(j,k)]]=beta[[paste(j,k)]] +
                        (Pst[i,j]-Pst[i,j]^2)*(t(Xij)%*%rij%*%t(rij)%*%Xij)/sigj^4 -
                            Pst[i,j]*(t(Xij)%*%Xij)/sigj^2
                    xi[[paste(j,k)]]=xi[[paste(j,k)]]+
                        (Pst[i,j]-Pst[i,j]^2)*(c(t(rij)%*%rij)/sigj^3-pi/sigj)^2+Pst[i,j]*(pi/sigj^2-3*c(t(rij)%*%rij)/sigj^4)
                    beta.xi[[paste(j,k)]] = beta.xi[[paste(j,k)]] +
                        ((Pst[i,j]-Pst[i,j]^2)*(c(t(rij)%*%rij)/sigj^3-pi/sigj)-Pst[i,j]*2/sigj)*t(Xij)%*%rij/sigj^2
                    if(X$G>1)
                    {
                        lambda[[paste(j,k)]]=lambda[[paste(j,k)]]+
                            Zi%*%t(Zi)*(Pst[i,j]*(1-2*P[i,j])+(2*P[i,j]^2-P[i,j])-Pst[i,j]^2 +2*Pst[i,j]*P[i,j] - P[i,j]^2 )
                        beta.lambda[[paste(j,k)]]=beta.lambda[[paste(j,k)]] +
                            (Pst[i,j]-Pst[i,j]^2)/sigj^2*t(Xij)%*%rij%*%t(Zi)
                        xi.lambda[[paste(j,k)]]=xi.lambda[[paste(j,k)]] +
                            (Pst[i,j]-Pst[i,j]^2)*(c(t(rij)%*%rij)/sigj^3 - pi/sigj)%*%t(Zi)

                    }


                }
                if(j<k)
                {
                    Xik=as.matrix(X$Xg[[k]][X$Yv[[i]],])
                    rik=c(X$Yl[[i]]- c(as.matrix(Xik)%*%X$param$mu[[k]]))
                    sigk=X$param$sigma[[k]][1]
                    beta[[paste(j,k)]]=beta[[paste(j,k)]] -
                        (Pst[i,j]*Pst[i,k])*(t(Xij)%*%rij%*%t(rik)%*%Xik)/sigj^2/sigk^2
                     xi[[paste(j,k)]]=xi[[paste(j,k)]] -
                         (Pst[i,j]*(c(t(rij)%*%rij)/sigj^3 -pi/sigj))*(Pst[i,k]*(c(t(rik)%*%rik)/sigk^3 -pi/sigk))
                    beta.xi[[paste(j,k)]] = beta.xi[[paste(j,k)]] +
                        -(Pst[i,j]*Pst[i,k])*(c(t(rik)%*%rik)/sigk^3-pi/sigk)*t(Xij)%*%rij/sigj^2
                    if(X$G>1)
                    {
                        lambda[[paste(j,k)]]=lambda[[paste(j,k)]]+
                            Zi%*%t(Zi)*( 2*P[i,j]*P[i,k]-(Pst[i,j]*P[i,k] + Pst[i,k]*P[i,j]) -(Pst[i,j]-P[i,j])*(Pst[i,k]-P[i,k]) )
                        beta.lambda[[paste(j,k)]]=beta.lambda[[paste(j,k)]] -
                             (Pst[i,j]*Pst[i,k])/sigj^2*t(Xij)%*%rij%*%t(Zi)
                        xi.lambda[[paste(j,k)]]=xi.lambda[[paste(j,k)]] -
                            (Pst[i,j]*Pst[i,k])*(c(t(rij)%*%rij)/sigj^3 - pi/sigj)%*%t(Zi)
                    }
                }
                if(i==X$N & j>k)
                {
                    beta[[paste(j,k)]]=t(beta[[paste(k,j)]])
                    xi[[paste(j,k)]]=t(xi[[paste(k,j)]])
                    beta.xi[[paste(j,k)]] = beta.xi[[paste(j,k)]] +
                        -(Pst[i,j]*Pst[i,k])*(c(t(rik)%*%rik)/sigk^3-pi/sigk)*t(Xij)%*%rij/sigj^2
                    if(X$G>1)
                    {
                        lambda[[paste(j,k)]]=t(lambda[[paste(k,j)]])
                        beta.lambda[[paste(j,k)]]=beta.lambda[[paste(j,k)]] -
                            (Pst[i,j]*Pst[i,k])/sigj^2*t(Xij)%*%rij%*%t(Zi)
                        xi.lambda[[paste(j,k)]]=xi.lambda[[paste(j,k)]] -
                            (Pst[i,j]*Pst[i,k])*(c(t(rij)%*%rij)/sigj^3 - pi/sigj)%*%t(Zi)
                    }
                }
            }
    }
    ### Construction de la matrice des BB'
    Hbb=numeric()
    for(i in 1:X$G)
    {
        line=numeric()
        for(j in 1:X$G)
        {
            line=cbind(line,beta[[paste(i,j)]])
        }
        Hbb=rbind(Hbb,line)
    }
    ### Construction de la matrice des xi.xi'
    Hxx=numeric()
    for(i in 1:X$G)
    {
        line=numeric()
        for(j in 1:X$G)
        {
            line=cbind(line,xi[[paste(i,j)]])
        }
        Hxx=rbind(Hxx,line)
    }
    ### Construction de la matrice des la la' si X$G>1
    if(X$G>1)
    {
        Hll=numeric()
        for(i in 2:X$G)
        {
            line=numeric()
            for(j in 2:X$G)
            {
                line=cbind(line,lambda[[paste(i,j)]])
            }
            Hll=rbind(Hll,line)
        }

    }
    ### Construction de la matrice des beta xi'
    Hbx=numeric()
    for(i in 1:X$G)
    {
        line=numeric()
        for(j in 1:X$G)
        {
            line=cbind(line,beta.xi[[paste(i,j)]])
        }
        Hbx=rbind(Hbx,line)
    }

    ## Construction de la matrice des beta lambda'
    if(X$G>1)
    {
        Hbl=numeric()
        for(i in 1:X$G)
        {
            line=numeric()
            for(j in 2:X$G)
            {
                line=cbind(line,beta.lambda[[paste(i,j)]])
            }
            Hbl=rbind(Hbl,line)
        }
    }

    ## Construction de la matrice des xi lambda'
    if(X$G>1)
    {
        Hxl=numeric()
        for(i in 1:X$G)
        {
            line=numeric()
            for(j in 2:X$G)
            {
                line=cbind(line,xi.lambda[[paste(i,j)]])
            }
            Hxl=rbind(Hxl,line)
        }
    }
    ### Construction de la matrice finale
    if(X$G==1)
    {
        H=rbind(cbind(Hbb,Hbx),cbind(t(Hbx),Hxx))
    }
    else if (X$G>1)
    {
        H=rbind(cbind(Hbb,Hbx,Hbl),cbind(t(Hbx),Hxx,Hxl),cbind(t(Hbl),t(Hxl),Hll))
    }
    -H
}


### calculons avec le gradient de la fonction de log[Likelihood]

IE.IND1=function(X)
{
    Hessian=matrix(0,X$pl+X$pm+X$pc,X$pl+X$pm+X$pc)
    Pst=post(X)
    if(X$G>1)
    {
        eta=X$Z%*%cbind(0,matrix(X$param$tau,ncol=X$G-1))
        P=exp(eta)/apply(exp(eta),1,sum)
    }
    for(i in 1:X$N)
    {
        beta=numeric()
        xi=numeric()
        la=numeric()
        for(j in 1:X$G)
        {
            Xij=as.matrix(X$Xg[[j]][X$Yv[[i]],])
            rij=c(X$Yl[[i]]- c(Xij%*%X$param$mu[[j]]))
            sig=X$param$sigma[[j]][1]
            beta=c(beta,Pst[i,j]*c(t(Xij)%*%rij)/sig^2)
            xi=c(xi,c(Pst[i,j]*c(t(rij)%*%rij/sig^3-X$pi[i])/sig))
            if(j>1)
            {
                la=c(la,X$Z[i,]*(Pst[i,j]-P[i,j]))
            }
        }
        Hessian=Hessian+c(beta,xi,la)%*%t(c(beta,xi,la))

    }
    Hessian
}
