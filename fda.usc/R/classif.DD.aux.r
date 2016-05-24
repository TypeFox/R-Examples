##############################
solve.ab <- function(a, b, tol = sqrt(.Machine$double.eps), LINPACK = FALSE, ...){
inv <- try(solve(a,b,tol, LINPACK,...),silent=TRUE)
if (class(inv)=="try-error") {   
  sv <-svd(a)
  inv <- drop((sv$v %*% diag(1/sv$d) %*% t(sv$u)) %*% b)
  #if (control$verbose) 
  #  warning("Inverse computed by SVD")
}
inv
}
##############################
quad.fit.opt <- function(index,dep,ind,tt){
  Df<-dep[,1]
  Dg<-dep[,2]
  A <- matrix(c(Df[index[1]],(Df[index[1]])^2, Df[index[2]],(Df[index[2]])^2),byrow=TRUE,2,2)
  w <- c(Dg[index[1]],Dg[index[2]])
  a <- solve.ab(A,w) 
  a02 <- a  
  if ( any(is.na(a02))) a <- c(1,1,1)
  mcr<-MCR0.p(a,dep,ind)
  a0 <- mcr0<- NULL 
  for (ii in tt) {
    a02 <- optim(a,AMCR0,dep=dep,ind=ind,tt=ii)$par
    mcr <- MCR0.p(a02,dep,ind)
    mcr0<-c(mcr0,mcr)
    a0<-rbind(a0,a02)
  }

  tt.opt <- tt[which.min(mcr0)]
  a0.opt <- a0[which.min(mcr0),]
#a   <- optim(a,AMCR0,dep=dep,ind=ind,tt=tt)$par  
#print(a)
#print(rbind(mcr,mcr2))
  return(list(min(mcr0), tt.opt, a0.opt))
}

##############################
cubic.fit.opt <- function(index,dep,ind,tt){
  Df<-dep[,1]
  Dg<-dep[,2]
  A <- matrix(c(Df[index[1]],(Df[index[1]])^2,(Df[index[1]])^3,
                Df[index[2]],(Df[index[2]])^2,(Df[index[2]])^3,
                Df[index[3]],(Df[index[3]])^2,(Df[index[3]])^3),byrow=TRUE,3,3)
  w <- c(Dg[index[1]],Dg[index[2]],Dg[index[3]])  
  a <- solve.ab(A,w)
  a02 <- a 
  if ( any(is.na(a02))) a <- c(1,1,1)
  mcr<-MCR0.p(a,dep,ind)
  a0 <- mcr0<- NULL 
  for (ii in tt) {
    a02 <- optim(a,AMCR0,dep=dep,ind=ind,tt=ii)$par
    mcr <- MCR0.p(a02,dep,ind)
    mcr0<-c(mcr0,mcr)
    a0<-rbind(a0,a02)
  }
  tt.opt <- tt[which.min(mcr0)]
  a0.opt <- a0[which.min(mcr0),]

  return(list(min(mcr0), tt.opt, a0.opt))
}

################################################################################
# Auxiliary functions for classif.DD
# File created by Manuel Oviedo de la Fuente  using code from paper:
# Li, J., P.C., Cuesta-Albertos, J.A. and Liu, R. 
# DD--Classifier: Nonparametric Classification Procedure Based on DD-plot. 
# Journal of the American Statistical Association (2012), Vol. 107, 737--753. 
##########################################################################################################################################
RR <- function(x,a){
   y <- 0
   kk <- length(a)
   for (i in 1:kk){
        y <- y+ a[i]*(x^i)
   }
   return(y)
}

AMCR0 <- function(a,dep,ind,tt){
  p<-sum(ind[,1])/nrow(ind)
  amcr <- p*mean(1/(1+exp(-tt*(dep[ind[,1],2]-sapply(dep[ind[,1],1],RR,a=a)))))+
        (1-p)*mean(1/(1+exp(-tt*(sapply(dep[ind[,2],1],RR,a=a)-dep[ind[,2],2]))))
  return(amcr)
}


quad.fit0 <- function(index,dep,ind){
   Df<-dep[,1]
   Dg<-dep[,2]
    A <- matrix(c(Df[index[1]],(Df[index[1]])^2, Df[index[2]],(Df[index[2]])^2),byrow=TRUE,2,2)
    w <- c(Dg[index[1]],Dg[index[2]])
    a <- try(solve(A,w),silent=TRUE)
    if (class(a)=="try-error") {
     sv<-svd(A)
     a<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u))%*%w)
#    warning("Inverse of sigma computed by SVD")
    }
    mcr<-MCR0.p(a,dep,ind)
    return(mcr)
}

cubic.fit0 <- function(index,dep,ind){
   Df<-dep[,1]
   Dg<-dep[,2]
    A <- matrix(c(Df[index[1]],(Df[index[1]])^2,(Df[index[1]])^3,
          Df[index[2]],(Df[index[2]])^2,(Df[index[2]])^3,
          Df[index[3]],(Df[index[3]])^2,(Df[index[3]])^3),byrow=TRUE,3,3)
    w <- c(Dg[index[1]],Dg[index[2]],Dg[index[3]])
    a <- try(solve(A,w),silent=TRUE)

    if (class(a)=="try-error") {
     sv<-svd(A)
     a<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u))%*%w)
 #    warning("Inverse of sigma computed by SVD")
    }
    mcr<-MCR0.p(a,dep,ind)
    return(mcr)
}

##########################################################################################################################################
# the function calculates the empirical misclassification rate from xx and yy based on the polynomial curve with coefficients being a in the DD-plot
# (the points on the curve are classified using KNN) NO

MCR0.p <- function(a,dep,ind){
 dm<-dim(ind)
 mis<-0
 p<-colMeans(ind)
 ahorro<-sapply(dep[ind[,1],1],RR,a=a)
 ahorro2<-sapply(dep[ind[,2],1],RR,a=a)
 for (i in 1:(dm[2]-1)) {
  mis<-p[i]*mean(ahorro<=dep[ind[,1],2])+(1-p[i])*mean(ahorro2>dep[ind[,2],2])
  }
 mis
 }


##########################################################################################################################################
# the function calculates the empirical misclassification rate from xx and yy based on the linear line with slope being k in the DD-plot
# (the points on the linear line are classified using KNN) NO

MCR0 <- function(k,dep,ind){
dm<-dim(ind)
mis<-0
p<-colMeans(ind)
for (i in 1:(dm[2]-1)) {
  mis<-mis+p[i]*mean((k*dep[ind[,i],1])<=dep[ind[,i],2])+(1-p[i])*mean((k*dep[-ind[,i],1])>dep[-ind[,i],2])
  }
 mis
 }
##########################################################################################################################################
##########################################################################################################################################
MCR.mof<-function(k,x,y,xx,yy,Dff,Dgg,n1,n2,nn1,nn2){
  p <- 0.5
  cls <- factor(c(rep(1,NROW(x)),rep(2,NROW(y))))  
  #print("MCR.mof")
  #print(length(cls))
  #print(dim(rbind(x,y)))
  xx0 <- as.matrix(xx[k*Dff[1:nn1]==Dgg[1:nn1],])
  yy0 <- as.matrix(yy[k*Dff[(nn1+1):(nn1+nn2)]==Dgg[(nn1+1):(nn1+nn2)],])
  if (ncol(xx0)==1) xx0 <- t(xx0)
  if (ncol(yy0)==1) yy0 <- t(yy0)
  mis.x <- 0
  mis.y <- 0
  
  if (nrow(xx0)!=0){
    
    #cls.x <- knn(rbind(x,y),xx0,cls,3)
    #mis.x <- sum(cls.x!=1)
  }
  if (nrow(yy0)!=0){
    #cls.y <- knn(rbind(x,y),yy0,cls,3)
    #mis.y <- sum(cls.y!=2)
  }
  mcr <- p*(sum(k*Dff[1:nn1]<Dgg[1:nn1])+mis.x)/nn1+
    (1-p)*(sum(k*Dff[(nn1+1):(nn1+nn2)]>Dgg[(nn1+1):(nn1+nn2)])+mis.y)/nn2
  return(mcr)
}

DD.cv.depth <- function(x,y,a0.2,a0.3,nam){
  
  n1 <- dim(x)[1]
  n2 <- dim(y)[1]  
  mis.cv.M <- mis.cv.D <- mis.cv.D3 <- mis.cv.D2 <- 0 
  for (l in 1:50){
    x.test <- x[1:(n1/50)+(l-1)*n1/50,]
    y.test <- y[1:(n2/50)+(l-1)*n2/50,]
    x.train <- x[-(1:(n1/50)+(l-1)*n1/50),]
    y.train <- y[-(1:(n2/50)+(l-1)*n2/50),]
    
    
    if (is.fdata(x)) {
      z.train <- c(x.train,y.train)
      z.test <- c(x.test,y.test)
      lst<-list(z.train,x.train)
      Df.train <- do.call(nam,lst)$dep
      lst<-list(z.train,y.train)
      Dg.train <- do.call(nam,lst)$dep
      lst<-list(z.test,x.train)
      Df.test <- do.call(nam,lst)$dep
      lst<-list(z.test,y.train)
      Dg.test <- do.call(nam,lst)$dep
      x.test<- x.test$data
      y.test<- y.test$data
      x.train<- x.train$data
      y.train<- y.train$data
      
    } else {
      z.train <- rbind(x.train,y.train)
      z.test <- rbind(x.test,y.test)
      lst<-list(z.train,x.train)
      Df.train <- do.call(nam,lst)$dep
      lst<-list(z.train,y.train)
      Dg.train <- do.call(nam,lst)$dep
      
      
      lst<-list(z.test,x.train)
      Df.test <- do.call(nam,lst)$dep
      lst<-list(z.test,y.train)
      Dg.test <- do.call(nam,lst)$dep
    }
    
    
    
    mis.cv.M  <- mis.cv.M + MCR.mof(1,x.train,y.train,x.test,y.test, Df.test,Dg.test,49/50*n1,49/50*n2,n1/50,n2/50)
    
    ##
    
    b <- Dg.train[Df.train!=0&Dg.train!=0]/Df.train[Df.train!=0&Dg.train!=0]
    b <- sort(b)
    m <- length(b)
    mis <- rep(0,m)
    
    mis <- sapply(b, MCR.mof,x=x.train,y=y.train,xx=x.train,yy=y.train,Dff=Df.train,Dgg=Dg.train,n1=49/50*n1,n2=49/50*n2,nn1=49/50*n1,nn2=49/50*n2)
    b0 <- min(b[which.min(mis)])
    
    mis.cv.D  <- mis.cv.D+ MCR.mof(b0,x.train,y.train,x.test,y.test,Df.test,Dg.test,49/50*n1,49/50*n2,n1/50,n2/50)
    
    ##
    dmax <- 1
    a2 <- optim(a0.2,AMCR,Dff=Df.train,Dgg=Dg.train,nn1=49/50*n1,nn2=49/50*n2,tt=100/dmax)$par
    mis.cv.D2 <- mis.cv.D2+MCR1.p(a2,x.train,y.train,x.test,y.test,Df.test,Dg.test,49/50*n1,49/50*n2,n1/50,n2/50)
    
    a2 <- optim(a0.3,AMCR,Dff=Df.train,Dgg=Dg.train,nn1=49/50*n1,nn2=49/50*n2,tt=100/dmax)$par
    mis.cv.D3 <- mis.cv.D3+MCR1.p(a2,x.train,y.train,x.test,y.test,Df.test,Dg.test,49/50*n1,49/50*n2,n1/50,n2/50)
  }
  return((c(mis.cv.M,mis.cv.D,mis.cv.D2,mis.cv.D3)+0.001*runif(4))/50)
}

AMCR <- function(a,Dff,Dgg,nn1,nn2,tt){
  p <- .5
  amcr <- p*sum(1/(1+exp(-tt*(Dgg[1:nn1]-sapply(Dff[1:nn1],RR,a=a)))))/nn1+
    (1-p)*sum(1/(1+exp(-tt*(sapply(Dff[(nn1+1):(nn1+nn2)],RR,a=a)-Dgg[(nn1+1):(nn1+nn2)]))))/nn2
  return(amcr)
}
##########################################################################################################################################     
# the function calculates the empirical misclassification rate from xx and yy based on the polynomial curve with coefficients being a in the DD-plot
# the points on the curve are classified using KNN 


MCR1.p <- function(a,x,y,xx,yy,Dff,Dgg,n1,n2,nn1,nn2){
  p <- .5
  cls <- factor(c(rep(1,n1),rep(2,n2)))
  xx0 <- as.matrix(xx[sapply(Dff[1:nn1],RR,a=a)==Dgg[1:nn1],])
  yy0 <- as.matrix(yy[sapply(Dff[(nn1+1):(nn1+nn2)],RR,a=a)==Dgg[(nn1+1):(nn1+nn2)],])
  if (ncol(xx0)==1) xx0 <- t(xx0)
  if (ncol(yy0)==1) yy0 <- t(yy0)
  mis.x <- 0
  mis.y <- 0
  if (nrow(xx0)!=0){
 #   cls.x <- knn(rbind(x,y),xx0,cls,3)
  #  mis.x <- sum(cls.x!=1)
  }
  if (nrow(yy0)!=0){
   # cls.y <- knn(rbind(x,y),yy0,cls,3)
   #  mis.y <- sum(cls.y!=2)
  }
  mcr.p <- p*(sum(sapply(Dff[1:nn1],RR,a=a)<Dgg[1:nn1])+mis.x)/nn1+
    (1-p)*(sum(sapply(Dff[(nn1+1):(nn1+nn2)],RR,a=a)>Dgg[(nn1+1):(nn1+nn2)])+mis.y)/nn2
  return(mcr.p)
}
