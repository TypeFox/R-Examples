"rda" <-
function(X,Y,scaling=1) {
# 
# Function RDA performs a redundancy analysis of the data in X
# and Y.
#
# scaling = 0 : use centred variables (X and Y)
# scaling = 1 : use centred and standardized variables (X and Y)
#
# Jan Graffelman
# Universitat Politecnica de Catalunya
# January 2004
#

n<-nrow(X)               # determine # of cases
p<-ncol(X)               # determine # of variables

nY<-nrow(Y)
q<-ncol(Y)

if (scaling==0) {
   Xa <- scale(X,scale=FALSE)
   Ya <- scale(Y,scale=FALSE)
}
else {
   if (scaling==1) {
      Xa<-scale(X)  
      Ya<-scale(Y)
   }
   else
      stop("rda: improper scaling parameter")
}

Rxx <- cor(X)
Ryy <- cor(Y)

B<-solve(t(Xa)%*%Xa)%*%t(Xa)%*%Ya

Yh<-Xa%*%B

pca.results <- princomp(Yh,cor=FALSE)
Fp <- pca.results$scores
Ds <- diag(sqrt(diag(var(Fp))))
Fs <- Fp%*%solve(Ds)
Gs <- pca.results$loadings
Gp <- Gs%*%Ds

la <- diag(var(Fp))
laf <- la/sum(la)
lac <- cumsum(laf)

decom <- rbind(la,laf,lac)

# It is important not to standardize Yh, and not to use cor=T. This will inflate the
# variance of Yh, and give different eigenvalues.

# Fs, Gp' biplot of fitted values

Fs <- Fs[,1:min(p,q)]
Fp <- Fp[,1:min(p,q)]
Gs <- Gs[,1:min(p,q)]
Gp <- Gp[,1:min(p,q)]

Gxs <- solve(t(Xa)%*%Xa)%*%t(Xa)%*%Fs
Gxp <- solve(t(Xa)%*%Xa)%*%t(Xa)%*%Fp

Gyp <- Gp
Gys <- Gs

# Gxs Gyp', Gxp Gys'  biplots of B (regression coefficients)

# goodness of fit of regression coefficients

decom <- decom[,1:min(p,q)]

# alternative computations doing SVD of B.

#W<-t(Xa)%*%Xa

#result <- svd(1/(sqrt(n-1))*half(W)%*%B)
#Gxxs <- sqrt(n-1)*mhalf(W)%*%result$u                    # = Gxs
#Gyyp <- result$v%*%diag(result$d)                        # = Gyp
#Gxxp <- sqrt(n-1)*mhalf(W)%*%result$u%*%diag(result$d)   # = Gxp
#Gyys <- result$v                                         # = Gys

#dd <- result$d*result$d
#dds <- cumsum(dd)
#ddf <- dd/sum(dd)
#ddc <- cumsum(ddf)

#decB <- rbind(dd,dds,ddf,ddc)

#res<-t(Gyyp)%*%Gyyp # =	D^2
#res<-t(Gyys)%*%Gyys # =	I
#res<-t(Gxxp)%*%Rxx%*%Gxxp # = D^2
#res<-t(Gxxs)%*%Rxx%*%Gxxs # = I

return(list(Yh=Yh,B=B,decom=decom,Fs=Fs,Gyp=Gyp,Fp=Fp,Gys=Gys,Gxs=Gxs,Gxp=Gxp))
}

