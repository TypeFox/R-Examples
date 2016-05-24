#################################################
### Fonctions reliee a l'estimation de modele ###
### avec une covariance de type CS.          ###
###                                         ###
### Analyste: Charles-Edouard Giguere      ###
#############################################

### On cree une liste de valeur qui pour une valeur de rho donnee retourne une liste de
### resultat intermediaire de facon a ne pas les recalculer par la suite.
###
### La fonction que l'on veut minimiser est la suivante :
### \sum_{i=1}^N[\tau_{ij}log|\Lambda_{ij}|] +  ( pfQ1 )
### N_jlog\sum_{i=1}^N[\tau_{ij](y_i-X_{ij}\hat{\beta}_j)'\Lambda_{ij}^{-1}(y_i-X_{ij}\hat{\beta}_j)]  ( pfQ2 )

pfQ.intermediate.CS1=function(rho,X,tau.ij,g)
{
    intermediate=list()

### Ce resultat est la premiere partie de la fonction profilee
### soit \sum_{i=1}^N\tau_{ij}log|\Lambda_{ij}|
    intermediate$pfQ1=sum(tau.ij*log((1-rho)^(X$pi-1)*(1+(X$pi-1)*rho)))

### Ce resultat est la derive de la premiere partie de la
### fonction profilee definie precedemment.
    intermediate$pfQ1.diff=-sum(tau.ij*(X$pi^2-X$pi)*rho/(1-rho)/(1+(X$pi-1)*rho))

### Ce resultat est la derive seconde de la premiere partie
### de la fonction profile definie precedemment.
    intermediate$pfQ1.hess=-sum(tau.ij*(X$pi^2-X$pi)*(1-(1-X$pi)*rho^2)/(1-rho)^2/(1+(X$pi-1)*rho)^2)

### Ce resultat genere les p inverses de la matrice Lambda pour les calculs subsequents.
    intermediate$Lambda.inv=list()
    for(i in 1:X$p)
        intermediate$Lambda.inv[[i]]=solve(matrix(rho,i,i)+(1-rho)*diag(i))

### Ce resultat genere l'estimation de beta pour cette valeur de rho
    dimXg=dim(X$Xg[[g]])[2]
    SXX=matrix(0,dimXg,dimXg)
    SXy=numeric(dimXg)
    for( i in 1:X$N)
    {
        SXX=SXX+tau.ij[i]*t(X$Xg[[g]][X$Yv[[i]],])%*%intermediate$Lambda.inv[[X$pi[i]]]%*%X$Xg[[g]][X$Yv[[i]],]
        SXy=SXy+c(tau.ij[i]*t(X$Xg[[g]][X$Yv[[i]],])%*%intermediate$Lambda.inv[[X$pi[i]]]%*%X$Yl[[i]])
    }
    intermediate$beta=c(solve(SXX,SXy))
### Ce resultat genere l'estimation de pfQ2
    S1=S2=S3=0;
    for(i in 1:X$N)
    {
        ri=c((X$Yl[[i]]-as.matrix(X$Xg[[g]][X$Yv[[i]],])%*%intermediate$beta))
        S1=S1 + c(tau.ij[i]*t(ri)%*%intermediate$Lambda.inv[[X$pi[i]]]%*%ri)
        S2=S2 + c(tau.ij[i]*t(ri)%*%intermediate$Lambda.inv[[X$pi[i]]]%*%
        (matrix(1,X$pi[i],X$pi[i])-diag(X$pi[i]))%*%
        intermediate$Lambda.inv[[X$pi[i]]]%*%ri)
        S3=S3 + 2*c(tau.ij[i]*t(ri)%*%intermediate$Lambda.inv[[X$pi[i]]]%*%
        (matrix(1,X$pi[i],X$pi[i])-diag(X$pi[i]))%*%
        intermediate$Lambda.inv[[X$pi[i]]]%*%
        (matrix(1,X$pi[i],X$pi[i])-diag(X$pi[i]))%*%
        intermediate$Lambda.inv[[X$pi[i]]]%*%ri)
    }
    intermediate$Ng=c(t(tau.ij)%*%X$pi)
    intermediate$sig2=S1/intermediate$Ng
    intermediate$pfQ2=intermediate$Ng*log(S1)
    intermediate$pfQ2.diff=-intermediate$Ng*S2/S1
    intermediate$pfQ2.hess=intermediate$Ng*(S3/S1-(S2/S1)^2)
    intermediate
}

estimloc.disp.CS1=function(X,rho,Post,iterlim,tol,iterEM)
{
    rho0=rho;
    rho1=numeric(X$G);
    sigma=beta=list()
    for(g in 1:X$G)
    {
        interm0=pfQ.intermediate.CS1(rho0[g],X,c(Post[,g]),g)
        for(iter in 1:iterlim)
        {
            if(iter==1)
                interm1=interm0
            else
                interm1=pfQ.intermediate.CS1(rho0[g],X,c(Post[,g]),g)
            rho1[g]=rho0[g]-(interm1$pfQ1.diff+interm1$pfQ2.diff)/(interm1$pfQ1.hess+interm1$pfQ2.hess)
            if(max(abs(rho1[g]-rho0[g]))<tol & rho1[g]>=0 & rho1[g]<1)
            {
                interm1=pfQ.intermediate.CS1(rho1[g],X,c(Post[,g]),g)
                sigma[[g]]=c(sqrt(interm1$sig2),rho1[g])
                beta[[g]]=c(interm1$beta)
                break
            }
            else if (rho1[g]>=1)
            {
                sigma[[g]]=c(sqrt(interm0$sig2),rho0[g])
                beta[[g]]=c(interm0$beta)
                warning(paste("At EM iteration:",iterEM,
                              "estimation of rho exceeded 1. Returned last good parameters."))
                break
            }
            else if (rho0[g]<0)
            {
                interm0=pfQ.intermediate.CS1(0,X,c(Post[,g]),g)
                sigma[[g]]=c(sqrt(interm0$sig2),0)
                beta[[g]]=c(interm0$beta)
                warning(paste("At EM iteration:",iterEM,
                              "estimation of rho was below 0. Rho is set at a value of 0."))
                break
            }
            else if(iter==iterlim)
            {
                sigma[[g]]=c(sqrt(interm1$sig2),rho1[g])
                beta[[g]]=c(interm1$beta)
                warning(paste("At EM iteration:",iterEM,
                              "estimation of covariance parameters did not reach convergence in group",g))
                break
            }
            else
            {
                rho0[g]=rho1[g]
                interm0=interm1;
            }
        }
    }
    return(list(beta,sigma))
}


##### Algorithme EM applique a des melanges multinormaux avec covariance CS et covariance inegale.


estimmmelnCS1=function(X,param,iterlim,tol)
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
        rho=numeric(X$G)
        for(g in 1:X$G)
        {
            rho[g]=sigma1[[g]][2]
        }
        loc.disp=estimloc.disp.CS1(X,rho,Post,iterlim,tol,iterationEM)
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


#### methode qui calcule la matrice hessienne des parametres du melange.
#### idealement ne rouler qu'une fois car cette methode peut s'averer un
#### peu lente.

I.CS1=function(X)
{
### on evalue la dimension de la matrice et on l'initialise.

    H=matrix(0,X$pl+X$pc+X$pm,X$pl+X$pc+X$pm)
    Pst=post(X)
### on evalue la sous-matrice des parametres de localisation.
    H.beta=list()
### on evalue la sous-matrice des parametres de covariance.
    H.sig=list()
    H.rho=list()
    H.rho.sig=list()
### on evalue la sous matrice des parametres de beta et covariance.
    H.beta.xi=list()

    if(X$G>1){
        H.lambda=list()
        eta=X$Z%*%matrix(X$param$tau,ncol=X$G-1)
        P=cbind(1,exp(eta))/apply(cbind(1,exp(eta)),1,sum)
        H.beta.lambda=list()
        H.xi.lambda=list()
    }
    for(i in 1:X$G)
    {
        for(j in i:X$G)
        {
            if(i==j) ### derive seconde de log L par theta_j^2
            {
                H.beta[[paste(j,j)]]=matrix(0,dim(X$Xg[[j]])[2],dim(X$Xg[[j]])[2])
                H.beta.xi[[paste(i,j)]]=matrix(0,dim(X$Xg[[j]])[2],2)
                H.sig[[paste(i,j)]]=0
                H.rho[[paste(i,j)]]=0
                H.rho.sig[[paste(i,j)]]=0
                SIG=cov.tsf(X$param$sigma[[j]],X$cov,X$p)
                Sig.invl=list()
                for(k in 1:X$p)
                {
                    Sig.invl[[k]]=solve(SIG[1:k,1:k])
                }
                if (X$G>1 & i>1)
                {
                    H.lambda[[paste(i,j)]]=0
                    H.beta.lambda[[paste(i,j)]]=matrix(0,dim(X$Xg[[j]])[2],dim(X$Z)[2])
                    H.xi.lambda[[paste(i,j)]]=matrix(0,2,dim(X$Z)[2])
                }
                for(k in 1:X$N)
                {
                    Xij=X$Xg[[j]][X$Yv[[k]],]
                    Sig=SIG[X$Yv[[k]],X$Yv[[k]]]
                    Sig.inv=Sig.invl[[X$pi[k]]]
                    rij=c(X$Yl[[k]]-as.matrix(Xij)%*%X$param$mu[[j]])
                    H.beta[[paste(j,j)]]= H.beta[[paste(j,j)]]+
                        Pst[k,j]*(1-Pst[k,j]) *
                            (t(Xij)%*%Sig.inv%*%rij%*%t(rij)%*%Sig.inv%*%Xij) -
                                Pst[k,j]*(t(Xij)%*%Sig.inv%*%Xij)

                    sig=X$param$sigma[[i]][1]
                    H.sig[[paste(i,j)]]=H.sig[[paste(i,j)]]+Pst[k,i]*
                        ((X$pi[k]^2+X$pi[k])*sig^(-2)-(2*X$pi[k]+3)*sig^(-2)*c(t(rij)%*%Sig.inv%*%rij)+
                         sig^(-2)*c((t(rij)%*%Sig.inv%*%rij)^2))-Pst[k,i]^2 *
                             (c(t(rij)%*%Sig.inv%*%rij)/sig-X$pi[k]/sig)^2
                    A=sig^2*(matrix(1,X$pi[k],X$pi[k])-diag(X$pi[k]))
                    H.rho[[paste(i,j)]]=H.rho[[paste(i,j)]] +
                        (Pst[k,i]-Pst[k,i]^2)*
                            (-1/2*sum(diag(Sig.inv%*%A)) + 1/2*c(t(rij)%*%Sig.inv%*%A%*%Sig.inv%*%rij))^2 +
                            Pst[k,i]*(-sum(diag(-Sig.inv%*%A%*%Sig.inv%*%A))/2 -
                                      c(t(rij)%*%Sig.inv%*%A%*%Sig.inv%*%A%*%Sig.inv%*%rij))
                    H.rho.sig[[paste(i,j)]]=H.rho.sig[[paste(i,j)]] +
                        (Pst[k,i]-Pst[k,i]^2)*(c(t(rij)%*%Sig.inv%*%rij)/sig-X$pi[k]/sig) *
                            (-1/2*sum(diag(Sig.inv%*%A))+1/2*c(t(rij)%*%Sig.inv%*%A%*%Sig.inv%*%rij)) -
                                Pst[k,i]/sig*(t(rij)%*%Sig.inv%*%A%*%Sig.inv%*%rij)
                    H.beta.xi[[paste(j,j)]]=H.beta.xi[[paste(j,j)]]+
                        cbind((Pst[k,i]*(c(t(rij)%*%Sig.inv%*%rij)/sig - X$pi[k]/sig -2/sig) -
                              Pst[k,i]^2*(c(t(rij)%*%Sig.inv%*%rij)/sig -X$pi[k]/sig))*t(Xij)%*%Sig.inv%*%rij ,
                              Pst[k,i]*(-2*t(Xij)%*%Sig.inv%*%A%*%Sig.inv%*%rij +
                                        (1-Pst[k,i])*(-sum(diag(Sig.inv%*%A))/2 +
                                                      c(t(rij)%*%Sig.inv%*%A%*%Sig.inv%*%rij)/2)*
                                        t(Xij)%*%Sig.inv%*%rij))
                    if(X$G>1 &i>1)
                    {
                        pi.ij=P[k,i]
                        H.lambda[[paste(i,j)]]=H.lambda[[paste(i,j)]] +
                            (X$Z[k,]%*%t(X$Z[k,]))*((Pst[k,i]*(1-2*pi.ij)+(2*pi.ij^2-pi.ij))-
                                                    (Pst[k,i]^2-2*pi.ij*Pst[k,i]+pi.ij^2))
                        H.beta.lambda[[paste(i,j)]]=H.beta.lambda[[paste(i,j)]] +
                            Pst[k,i]*(1-pi.ij)*t(Xij)%*%Sig.inv%*%rij%*%t(X$Z[k,]) - (Pst[k,i]^2-Pst[k,i]*pi.ij)*
                                t(Xij)%*%Sig.inv%*%rij%*%t(X$Z[k,])
                        H.xi.lambda[[paste(i,j)]]=H.xi.lambda[[paste(i,j)]] +
                            rbind(
                                  (Pst[k,i]*(1-pi.ij)-Pst[k,i]*(Pst[k,j]-pi.ij))*(c(t(rij)%*%Sig.inv%*%rij)/
                                                                                  sig-X$pi[k]/sig)*t(X$Z[k,])
                                  ,
                                  (Pst[k,i]*(1-pi.ij)-Pst[k,i]*(Pst[k,j]-pi.ij))*
                                (-1/2*sum(diag(Sig.inv%*%A))+1/2*c(t(rij)%*%Sig.inv%*%A%*%Sig.inv%*%rij))*
                                t(X$Z[k,])
                                  )

                    }
                }
            }
            else if(i!=j) ### derive seconde de log L par theta_i x theta_j
            {
                H.beta[[paste(i,j)]]=matrix(0,dim(X$Xg[[i]])[2],dim(X$Xg[[j]])[2])
                H.beta.xi[[paste(i,j)]]=0
                H.beta.xi[[paste(j,i)]]=0
                H.sig[[paste(i,j)]]=0
                H.rho[[paste(i,j)]]=0
                H.rho.sig[[paste(i,j)]]=0
                H.rho.sig[[paste(j,i)]]=0
                if(X$G>1 & i>1)
                    H.lambda[[paste(i,j)]]=0
                if(X$G>1 & j>1)
                {
                    H.beta.lambda[[paste(i,j)]]=0
                    H.xi.lambda[[paste(i,j)]]=matrix(0,2,dim(X$Z)[2])

                }
                if(X$G>1 & i>1)
                {
                    H.beta.lambda[[paste(j,i)]]=0
                    H.xi.lambda[[paste(j,i)]]=matrix(0,2,dim(X$Z)[2])
                }
                SIGi=cov.tsf(X$param$sigma[[i]],X$cov,X$p)
                SIGj=cov.tsf(X$param$sigma[[j]],X$cov,X$p)

                for(k in 1:X$N)
                {
                    Xij.i=X$Xg[[i]][X$Yv[[k]],]
                    Xij.j=X$Xg[[j]][X$Yv[[k]],]
                    Sig.inv.i=solve(SIGi[X$Yv[[k]],X$Yv[[k]]])
                    Sig.inv.j=solve(SIGj[X$Yv[[k]],X$Yv[[k]]])
                    rij.i=c(X$Yl[[k]]-as.matrix(Xij.i)%*%X$param$mu[[i]])
                    rij.j=c(X$Yl[[k]]-as.matrix(Xij.j)%*%X$param$mu[[j]])
                    sig.i=X$param$sigma[[i]][1]
                    sig.j=X$param$sigma[[j]][1]
                    A.i=sig.i^2*(matrix(1,X$pi[k],X$pi[k])-diag(X$pi[k]))
                    A.j=sig.j^2*(matrix(1,X$pi[k],X$pi[k])-diag(X$pi[k]))
                    H.beta[[paste(i,j)]]=H.beta[[paste(i,j)]]-
                        Pst[k,i]*Pst[k,j]*(t(Xij.i)%*%Sig.inv.i%*%rij.i%*%t(rij.j)%*%Sig.inv.j%*%Xij.j)
                    H.sig[[paste(i,j)]]=H.sig[[paste(i,j)]]-
                        Pst[k,i]*Pst[k,j]*(c(t(rij.i)%*%Sig.inv.i%*%rij.i)/sig.i-X$pi[k]/sig.i)*
                            (c(t(rij.j)%*%Sig.inv.j%*%rij.j)/sig.j-X$pi[k]/sig.j)
                    H.rho[[paste(i,j)]]=H.rho[[paste(i,j)]] -
                        Pst[k,i]*Pst[k,j]*(-1/2*sum(diag(Sig.inv.i%*%A.i)) +
                                           1/2*c(t(rij.i)%*%Sig.inv.i%*%A.i%*%Sig.inv.i%*%rij.i)) *
                            (-1/2*sum(diag(Sig.inv.j%*%A.j)) +
                             1/2*c(t(rij.j)%*%Sig.inv.j%*%A.j%*%Sig.inv.j%*%rij.j))
                    H.rho.sig[[paste(i,j)]]=H.rho.sig[[paste(i,j)]] -
                        Pst[k,i]*Pst[k,j]*(c(t(rij.i)%*%Sig.inv.i%*%rij.i)/sig.i-X$pi[k]/sig.i) *
                            (-1/2*sum(diag(Sig.inv.j%*%A.j)) +
                             1/2*c(t(rij.j)%*%Sig.inv.j%*%A.j%*%Sig.inv.j%*%rij.j))
                    H.rho.sig[[paste(j,i)]]=H.rho.sig[[paste(j,i)]] -
                        Pst[k,i]*Pst[k,j]*(c(t(rij.j)%*%Sig.inv.j%*%rij.j)/sig.j-X$pi[k]/sig.j) *
                            (-1/2*sum(diag(Sig.inv.i%*%A.i)) +
                             1/2*c(t(rij.i)%*%Sig.inv.i%*%A.i%*%Sig.inv.i%*%rij.i))
                    pi.i=P[k,i]
                    pi.j=P[k,j]
                    if(X$G>1 & i>1)
                        H.lambda[[paste(i,j)]]=H.lambda[[paste(i,j)]] +
                            (X$Z[k,]%*%t(X$Z[k,]))*((2*pi.i*pi.j-(Pst[k,i]*pi.j+Pst[k,j]*pi.i)) -
                                                    (Pst[k,i]-pi.i)*(Pst[k,j]-pi.j))
                    H.beta.xi[[paste(i,j)]]=H.beta.xi[[paste(i,j)]]+
                        cbind(-Pst[k,i]*Pst[k,j]*(c(t(rij.j)%*%Sig.inv.j%*%rij.j)/sig.j -
                                                  X$pi[k]/sig.j)*t(Xij.i)%*%Sig.inv.i%*%rij.i ,
                              -Pst[k,i]*Pst[k,j]*
                              (-sum(diag(Sig.inv.j%*%A.j))/2 +
                               c(t(rij.j)%*%Sig.inv.j%*%A.j%*%Sig.inv.j%*%rij.j)/2)*
                              t(Xij.i)%*%Sig.inv.i%*%rij.i)
                     H.beta.xi[[paste(j,i)]]=H.beta.xi[[paste(j,i)]]+
                        cbind(-Pst[k,i]*Pst[k,j]*(c(t(rij.i)%*%Sig.inv.i%*%rij.i)/sig.i -
                                                  X$pi[k]/sig.i)*t(Xij.j)%*%Sig.inv.j%*%rij.j ,
                              -Pst[k,i]*Pst[k,j]*(-sum(diag(Sig.inv.i%*%A.i))/2 +
                                                  c(t(rij.i)%*%Sig.inv.i%*%A.i%*%Sig.inv.i%*%rij.i)/2)*
                              t(Xij.j)%*%Sig.inv.j%*%rij.j)
                    if(X$G>1 & j>1)
                    {
                        H.beta.lambda[[paste(i,j)]]=H.beta.lambda[[paste(i,j)]] +
                            -Pst[k,i]*pi.j*t(Xij.i)%*%Sig.inv.i%*%rij.i%*%t(X$Z[k,]) -
                                Pst[k,i]*(Pst[k,j]-pi.j)*t(Xij.i)%*%Sig.inv.i%*%rij.i%*%t(X$Z[k,])
                        H.xi.lambda[[paste(i,j)]]=H.xi.lambda[[paste(i,j)]] +
                            rbind(
                                  (-Pst[k,i]*pi.j-Pst[k,i]*(Pst[k,j]-pi.j))*(c(t(rij.i)%*%Sig.inv.i%*%rij.i)/sig.i
                                                                             - X$pi[k]/sig.i)*t(X$Z[k,]),
                                  (-Pst[k,i]*pi.j-Pst[k,i]*(Pst[k,j]-pi.j))*
                                  (-1/2*sum(diag(Sig.inv.i%*%A.i))+
                                   1/2*c(t(rij.i)%*%Sig.inv.i%*%A.i%*%Sig.inv.i%*%rij.i))*t(X$Z[k,])
                                  )

                    }
                    if(X$G>1 & i>1)
                    {
                        H.beta.lambda[[paste(j,i)]]=H.beta.lambda[[paste(j,i)]] +
                            -Pst[k,j]*pi.i*t(Xij.j)%*%Sig.inv.j%*%rij.j%*%t(X$Z[k,]) -
                                Pst[k,j]*(Pst[k,i]-pi.i)*t(Xij.j)%*%Sig.inv.j%*%rij.j%*%t(X$Z[k,])
                        H.xi.lambda[[paste(j,i)]]=H.xi.lambda[[paste(j,i)]] +
                            rbind(
                                  (-Pst[k,j]*pi.i-Pst[k,j]*(Pst[k,i]-pi.i))*
                                (c(t(rij.j)%*%Sig.inv.j%*%rij.j)/sig.j -
                                 X$pi[k]/sig.j)*t(X$Z[k,])
                                  ,
                                  (-Pst[k,j]*pi.i-Pst[k,j]*(Pst[k,i]-pi.i))*
                                  (-1/2*sum(diag(Sig.inv.j%*%A.j)) +
                                   1/2*c(t(rij.j)%*%Sig.inv.j%*%A.j%*%Sig.inv.j%*%rij.j))*t(X$Z[k,])
                                  )

                    }
                }
                if(X$G>1 & i>1)
                    H.lambda[[paste(j,i)]] = t(H.lambda[[paste(i,j)]])
                H.beta[[paste(j,i)]] = t(H.beta[[paste(i,j)]])

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
            line=cbind(line,H.beta[[paste(i,j)]])
        }
        Hbb=rbind(Hbb,line)
    }
### Construction de la matrice des xi.xi'

    Hxx=matrix(0,X$pc,X$pc)
    for(gl in 1:X$G) ### groupe dans la ligne
    {
        for( vl in 1:2)  ### 1=sig 2=rho identification de la variable de la ligne
        {
            for( gc in 1:X$G) ### groupe dans la colonne
            {
                for( vc in 1:2) ### 1=sig 2=rho identification de la variable de la colonne.
                {
                    if( vl==vc & vl==1)
                    {
                        if(gl<=gc)
                            Hxx[2*(gl-1)+vl,2*(gc-1)+vc]=H.sig[[paste(gl,gc)]]
                        else
                            Hxx[2*(gl-1)+vl,2*(gc-1)+vc]=H.sig[[paste(gc,gl)]]
                    }
                    else if( vl==vc & vl==2)
                    {
                        if(gl<=gc)
                            Hxx[2*(gl-1)+vl,2*(gc-1)+vc]=H.rho[[paste(gl,gc)]]
                        else
                            Hxx[2*(gl-1)+vl,2*(gc-1)+vc]=H.rho[[paste(gc,gl)]]

                    }
                    else if (vl!=vc)
                    {
                        if(vl<vc)
                            Hxx[2*(gl-1)+vl,2*(gc-1)+vc]=H.rho.sig[[paste(gl,gc)]]
                        else
                            Hxx[2*(gl-1)+vl,2*(gc-1)+vc]=H.rho.sig[[paste(gc,gl)]]
                    }
                }
            }
        }
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
                line=cbind(line,H.lambda[[paste(i,j)]])
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
            line=cbind(line,H.beta.xi[[paste(i,j)]])
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
                line=cbind(line,H.beta.lambda[[paste(i,j)]])
            }
            Hbl=rbind(Hbl,line)
        }
    }

    ## Construction de la matrice des beta lambda'
    if(X$G>1)
    {
        Hxl=numeric()
        for(i in 1:X$G)
        {
            line=numeric()
            for(j in 2:X$G)
            {
                line=cbind(line,H.xi.lambda[[paste(i,j)]])
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

IE.CS1=function(X)
{
    Hessian=matrix(0,X$pl+X$pm+X$pc,X$pl+X$pm+X$pc)
    Pst=post(X)
    if(X$G>1)
    {
        eta=X$Z%*%cbind(0,matrix(X$param$tau,ncol=X$G-1))
        P=exp(eta)/apply(exp(eta),1,sum)
    }
    SIG=Sig.invl=list()
    for(g in 1:X$G)
    {
        SIG[[g]]=cov.tsf(X$param$sigma[[g]],X$cov,X$p)
        for(k in 1:X$p)
        {
            Sig.invl[[paste(g,k)]]=solve(SIG[[g]][1:k,1:k])
        }
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
            Sig.inv=Sig.invl[[paste(j,X$pi[i])]]
            A=matrix(1,X$pi[i],X$pi[i])-diag(X$pi[i])
            beta=c(beta,Pst[i,j]*c(t(Xij)%*%Sig.inv%*%rij))
            xi=c(xi,c(Pst[i,j]*c(t(rij)%*%Sig.inv%*%rij-X$pi[i])/sig),
            c(Pst[i,j]*sig^2/2*(c(t(rij)%*%Sig.inv%*%A%*%Sig.inv%*%rij)-sum(diag(Sig.inv%*%A)))))
            if(j>1)
            {
                la=c(la,X$Z[i,]*(Pst[i,j]-P[i,j]))
            }
        }
        Hessian=Hessian+c(beta,xi,la)%*%t(c(beta,xi,la))

    }
    Hessian
}
