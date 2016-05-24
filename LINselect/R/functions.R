############################################################################

ProjY <- function(y,V,d,min.vp=10**(-8)) {
  # ---------------------------------------------------------------
  # COMPUTE THE PROJECTION OF y OVER THE SPACE SPANNED BY V
  # called by VARselect
  # INPUT
  # V: n x d matrix
  # y: n dim. vector
  # d: integer
  # min.vp: (small) poisitive real number
  # OUTPUT
  # Proj: n dim. vector
  # CALL:
  # ---------------------------------------------------------------
  # SH : 19/02/2013 modifs
  M <- t(V)%*%V
  svdM <- svd(M,nu=0,nv=0)
  s1 <- svdM$d[1]
  if (s1==0) {
    Proj <- y*0
    rgV <- 0
  }
  if (s1!=0) {
    Mn <- M/s1
    svdMn <- svd(Mn,nu=0,nv=0)
    rgV <- sum(svdMn$d > min.vp)
    if (rgV == 0) {
      Proj <- y*0
    }
    if (rgV == d) {
      Proj <- V%*%solve(M,t(V)%*%y)
    }
    if (rgV==1) {
      Proj <-1/s1*V%*%t(V)%*%y
    }
    if ((rgV > 1)&(rgV < d)) {
      svdMn <-  svd(Mn)
      invM <-
        svdMn$v[,1:rgV]%*%diag(1/svdMn$d[1:rgV])%*%t(svdMn$v[,1:rgV])/s1
      Proj <- V%*%invM%*%t(V)%*%y
    }
  }
  return(list(Proj=Proj,rg=rgV))
} #  fin ProjY

# ---------------------------------------------------
Proj <- function(V,d,min.vp=10**(-8)) {
  # ---------------------------------------------------------------
  # COMPUTE THE PROJECTION matrix OVER THE SPACE SPANNED BY V
  # MODIFIE PAR JEHANNE
  # INPUT
  # V: n x d matrix
  # d: integer
  # min.vp: (small) positive real number
  # OUTPUT
  # Proj: n x n  dim. matrix
  # CALL:
  # ---------------------------------------------------------------
  # SH : 19/02/2013 ajout de Proj
  M <- t(V)%*%V
  svdM <- svd(M,nu=0,nv=0)
  s1<-svdM$d[1]
  if (s1==0) {
    Proj<-matrix(0,ncol=0,nrow=0)
    rgV<-0
  }
  if (s1!=0) {
    Mn<-M/s1
    svdMn<-svd(Mn,nu=0,nv=0)
    rgV <- sum(svdMn$d > min.vp)
    if (rgV == d) {
      Proj <- V%*%solve(M)%*%t(V)
    }
    if (rgV==1) {
      Proj <-1/s1*V%*%t(V)
    }
    if ((rgV>1)&(rgV<d)) {
      svdMn <- svd(Mn)
      invM <-svdMn$v[,1:rgV]%*%diag(1/svdMn$d[1:rgV])%*%t(svdMn$v[,1:rgV])/s1
      Proj <- V%*%invM%*%t(V)
    }
    if (rgV==0) {
      Proj <- matrix(0,ncol=0,nrow=0)
    }
  }
  return(list(Proj=Proj,rg=rgV))
} #  fin Proj
 


# ---------------------------------------------------
active<-function(res,in_s) {
 # Calcule l'active set le long de la solution res (rendu d'enet) et
 #         en renumerotant selon in_s (vecteur) 
 # rend une liste de vecteurs, active set a  chaque etape de la
 #     trajectoire (indice des variables selectionnees) 
p<-length(res$penalty)
act<-list(NULL)
for (i in 1:p) {
	ind<-which(res$beta.pure[i,]!=0)
	act[[i]]<-in_s[res$allset[ind]]
}
return(act)
} # fin active

# ---------------------------------------------------
estime2<-function(Y,X,It,Il,act) {
 # Y est un vecteur de taille n (observations)
 # X est une matrice de taille n*p (covariables)
 # It est un vecteur contenant des indices entiers appartenant a
 #          1,..,n (indices de l'echantillon test)
 # Il est un vecteur contenant des indices entiers appartenant a
 #          1,..,n (indices de l'echantillon d'apprentissage)
  
 # act est une liste de vecteurs contenant des indices entiers
 #      appartenant a 1,..,p (indices des variables considerees comme
 #      explicatives, active set) 
 # estime2 rend est, la liste de l'erreur de prediction de Y[I]
  T<-length(act)
  est<-list(NULL)
for (i in 1:T) {
  if (length(act[[i]])>0) {
    beta <- lm(Y[Il]~X[Il,act[[i]]]-1)$coef
    if (length(beta)==1) Ypred <- X[It,act[[i]]]*beta
    if (length(beta)>1)  Ypred <- X[It,act[[i]]]%*%beta
  }
  else {
    Ypred<-rep(0,length(It))
  }
  est[[i]]<-sum((Y[It]-Ypred)^2)/length(It)
}
return(est)
}

# ---------------------------------------------------
indice<-function(l,lambda) {
 # si lambda est un vecteur (suppose decroissant) donne l'indice du
 #      plus petit element de lambda superieur ou egal a   l (nombre) 
  compare<-which(lambda>=l)
  if (length(compare)>0) {
    ind<-compare[length(compare)]
  }
  else {
    ind<-1
  }
  return(ind)
} # fin indice

# ---------------------------------------------------
enleve_var_0<-function(X,I) {
 # X matrice de taille n*p
 # I vecteur inclus dans 1:n
 # rend in_s un vecteur inclus dans 1:p contient le numero des
 #         variables explicatives de variance non nulle sur les observations I 
  s<-scale(X[I,])
  in_s<-which(attr(s,"scaled:scale")!=0)
  return(in_s)
} # fin enleve_var_0

#--------------------------
calc.fhat <- function(Y,X,BetaHat,Nmod,n,p) {
   # compute the projection of Y onto the spaces of selected variables
   # called by VARselect
   # Y n vector
   # X : n X p matrix
   # BetaHat : matrix with Nmod rows
   # Nmod, n, p integers
   # return the projection of Y onto the spaces of selected variables
   # and the rank of the spaces
   f.proj <- matrix(0,nrow=n,ncol=Nmod)
   r.proj <- rep(0, Nmod)
   m.lasso <- list(NULL)
   for (h in 1:Nmod) {
     mod <- (1:p)[BetaHat[h,]!=0]
     m.lasso[[h]] <- (1:p)[BetaHat[h,]!=0]
     if (length(mod)==0) {
       f.proj[,h] <- rep(mean(Y),n)
       r.proj[h] <- 1
       m.lasso[[h]] <- "Intercept"
     }
     if (length(mod)!=0) {
       Pr <- ProjY(Y,cbind(rep(1,n),X[,mod]),length(mod)+1)
       f.proj[,h] <- Pr$Proj
       r.proj[h] <- Pr$rg
     }
   }
   return(list(fhat=f.proj,r=r.proj,m=m.lasso))
 } # fin calc.fhat
 
EstSelect5 <- function(Y, X, fHat, IHat, rg, pen) {
  # Linear regression model : f=X beta
  # projection estimator on S_m
  # The approximate space = S_m
  #
  # INPUT
  # Y observations. vector n
  # X covariates. matrix nxp
  # fHat : matrix nxdimH.
  #        fHat[,h] projection of Y on a linear space IHat[[h]]+cste with
  #        dimension rg[h]
  # IHat : list 
  # rg : vector dimH.
  # pen : vector of length Dmax+1
  # 
  # OUTPUT : list with components
  # crit : vector with dimension dimH
  #        crit[h] : value of the criteria for IHat[,indS[l]]
  # lHat : integer. indice of the chosen estimator
  # fHat : vector of length n. estimated value of fHat[,lHat]
  #   
  n <- length(Y)
  dimH <- dim(fHat)[2]
  crit <- rep(0,dimH)
  for (h in 1:dimH) {
    N1 <- sum((Y-fHat[,h])**2)
    # pen(1) correspond a D=0
    SCRpen <- pen[rg[h]+1]*N1/(n-rg[h])
    crit[h] <- N1 + SCRpen
  }
  lHat <- which.min(crit)
  return(list(crit=crit,
              lHat=lHat,fHat=fHat[,lHat],
               mHat=IHat[[lHat]]))
}#  fin EstSelect5
