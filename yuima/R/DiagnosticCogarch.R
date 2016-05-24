# We write a Diagnostic function that evaluates the following quantity
Diagnostic.Cogarch <- function(yuima.cogarch, param = list(), 
                               matrixS = NULL , mu = 1,
                               display =TRUE){
  
  
  if(missing(yuima.cogarch))
    yuima.stop("yuima.cogarch or yuima object is missing.")
  
  if((length(param)==0 && !is(yuima.cogarch, "cogarch.gmm"))){
      yuima.stop("missing values parameters")
  }
  
  
  if(is(yuima.cogarch,"yuima")){
    model<-yuima.cogarch@model
  }else{
    if(is(yuima.cogarch,"yuima.cogarch")){
      model<-yuima.cogarch
    }else{
      if(is(yuima.cogarch,"cogarch.gmm")){
        model<-yuima.cogarch@model
        if(length(param)==0){
          param<-coef(yuima.cogarch)
        }
      }
    }
  }
  
  if(!is.COGARCH(model)){
    yuima.warn("The model does not belong to the class yuima.cogarch")
  }
  
  info <- model@info
  numb.ar <- info@q
  ar.name <- paste(info@ar.par,c(numb.ar:1),sep="")
  numb.ma <- info@p
  ma.name <- paste(info@ma.par,c(1:numb.ma),sep="")
  loc.par <- info@loc.par
  
  if(is.list(param)){
    param <- unlist(param)
  }  
  
  xinit.name0 <- model@parameter@xinit
  idx <- na.omit(match(c(loc.par, ma.name), xinit.name0))
  xinit.name <- xinit.name0[-idx]

  
  fullcoeff <- c(ar.name, ma.name, loc.par,xinit.name)
  nm<-names(param)
  oo <- match(nm, fullcoeff)
  if(!is(yuima.cogarch,"cogarch.gmm")){
    if(length(na.omit(oo))!=length(fullcoeff))
      yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima.cogarch model")
  }
  acoeff <- param[ma.name]
  b <- param[ar.name]
  cost<- param[loc.par]
  # Check Strictly Stationary condition
  Amatr<-MatrixA(b[c(info@q:1)])
  if(is.null(matrixS)){
    Dum <- yuima.carma.eigen(Amatr)
    matrixS <- Dum$vectors
    lambda.eig <- Dum$values
  }else{
    if(!is.matrix(matrixS)){
      if(dim(matrixS)[1]!=dim(matrixS)[2]){
        yuima.stop("matrixS must be an object of class matrix")
      }else{
        if(dim(matrixS)[1]!=numb.ar)
          yuima.stop("matrixS must be an square matrix. Its row number is equal to the dimension of autoregressive coefficients")
      }
    }
    lambda.eig<-diag(solve(matrixS)%*%Amatr%*%matrixS)
  }
  lambda1<- max(Re(lambda.eig)) # we find lambda1
  ev.dum<-matrix(0,1,info@q)
  ev.dum[1,info@q] <- 1
  av.dum <-matrix(0,info@q,1)
  av.dum[c(1:info@p),1]<-acoeff 
  if(display==TRUE){
    cat(paste0("\n COGARCH(",info@p,info@q,") model \n"))
  }
  matrforck<-solve(matrixS)%*%(av.dum%*%ev.dum)%*%matrixS
  if(is.complex(matrforck)){
    A2Matrix<-Conj(t(matrforck))%*%matrforck
    SpectNorm<-max(sqrt(as.numeric(abs(eigen(A2Matrix)$values))))
  }else{
    A2Matrix<-t(matrforck)%*%matrforck
    SpectNorm<-max(sqrt(abs(eigen(A2Matrix)$values)))
  }
  if(SpectNorm*mu < -lambda1){
    if(display==TRUE){
      cat("\n The process is strictly stationary\n The unconditional first moment of the Variance process exists \n")
    }
    res.stationarity <- TRUE
  }else{
    if(display==TRUE){
      cat("\n We are not able to establish if the process is stationary \n
          Try a new S matrix such that the Companion matrix is diagonalizable")
    }
    res.stationarity <- "Try a different S matrix"
  }
  # Check the positivity
  massage <- "\n We are not able to establish the positivity of the Variance \n"
  res.pos <- " "
  if(info@p==1){
    if(is.numeric(lambda.eig) && all(lambda.eig<0)){
      massage <- "\n the Variance is a positive process \n"
      res.pos <- TRUE
    }
    if(is.complex(lambda.eig) && all(Im(lambda.eig)==0) && all(Re(lambda.eig) < 0)){
      massage <- "\n the Variance is a positive process \n"
      res.pos <- TRUE
    } 
  }
  if((info@q==2 && info@p==2 )){
    if(is.numeric(lambda.eig)|| (is.complex(lambda.eig) && all(Im(lambda.eig)==0))){
      if(acoeff[1]>=-acoeff[2]*lambda1){
        massage <- "\n the Variance is a positive process \n"
        res.pos <- TRUE
      }else{
        massage <- "\n the Variance process is not strictly positive \n"
        res.pos <- FALSE
      }
    }
  }
  if(info@p>=2 && info@q>info@p){
    gamma.eig<-polyroot(acoeff)
    if(all(Im(gamma.eig)==0)){
      gamm.eig<-as.numeric(gamma.eig)
    }
    if(is.numeric(lambda.eig) && all(is.numeric(gamma.eig)<0) && all(is.numeric(lambda.eig)<0)){
      if(cumsum(sort(lambda.eig[c(1:(info@p-1))]))>=cumsum(sort(gamma.eig[c(1:(info@p-1))]))){
        massage <- "\n The Variance process is strictly positive. \n"
      }
    }
  }
  if(display==TRUE){
    cat(massage)
  }
  
  res1<-StationaryMoments(cost,b,acoeff,mu)
  res <- list(meanVarianceProc=res1$ExpVar,
              meanStateVariable=res1$ExpStatVar,
              stationary = res.stationarity,
              positivity = res.pos)
  return(res)

}

yuima.norm2<-function(x){
    res<-sqrt(max(eigen(t(x) %*% x)$values))
    return(res)
}

StationaryMoments<-function(cost,b,acoeff,mu=1){
  # We obtain stationary mean of State process
  # stationary mean of the variance 
  # E\left(Y\right)=-a_{0}m_{2}\left(A+m_{2}ea'\right)^{-1}e
  q<-length(b) 
  a <- e <- matrix(0,nrow=q,ncol=1)
  e[q,1] <- 1
  p<-length(acoeff)
  a[1:p,1] <- acoeff
  B_tilde <- MatrixA(b[c(q:1)])+mu*e%*%t(a)
  
  if(q>1){
    invB<-rbind(c(-B_tilde[q,-1],1)/B_tilde[q,1],cbind(diag(q-1),matrix(0,q-1,1)))
  }else{invB<-1/B_tilde}
  
  ExpStatVar <- -cost*mu*invB%*%e
  ExpVar <- cost+t(a)%*%ExpStatVar
  res <- list(ExpVar=ExpVar, ExpStatVar=ExpStatVar)
  return(res)
}

