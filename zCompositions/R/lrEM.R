lrEM <- function(X,label=NULL,dl=NULL,rob=FALSE,ini.cov=c("complete.obs","multRepl"),delta=0.65,tolerance=0.0001,
         max.iter=50,rlm.maxit=150,suppress.print=FALSE){
  
  if (is.character(dl)) stop("dl must be a numeric vector or matrix")
  if (is.vector(dl)) dl <- matrix(dl,nrow=1)
  
  if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a data matrix")
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)){
    if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    if (label!=0 & any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
    if (any(is.na(X))) stop(paste("NA values not labelled as censored values were found in the data set"))
  }
  if (is.na(label)){
    if (any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
    if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
  }
  if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
  if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
  
  ini.cov <- match.arg(ini.cov)
  
  lm.sweep <- function(M,C,varobs){
    
    sweep.matrix <- function(A,ind){ 
      
      nn <- nrow(A); D <- ncol(A)
      S <- A
      
      for (j in ind){                       
        S[j,j] <- -1/A[j,j]     
        for (i in 1:D) {           
          if (i != j){
            S[i,j] <- -A[i,j]*S[j,j]
            S[j,i] <- S[i,j]
          }
        }
        for (i in 1:D){           
          if (i != j){
            for (k  in 1:D){
              if (k != j){
                S[i,k] <- A[i,k] - S[i,j]*A[j,k]
                S[k,i] <- S[i,k]
              }
            }
          }
        }
        A <- S 
      }
      return(A)
    }  
    
    D <- length(M)                      
    q <- length(varobs)                     
    i <- rep(1,D)                  
    i[varobs] <- i[varobs]-1
    dep <- which(i!=0)
    ndep <- length(dep)                 
    
    A <- matrix(0,D+1,D+1)             
    A[1,1] <- -1                    
    A[1,2:(D+1)] <- M                
    A[2:(D+1),1] <- matrix(M,ncol=1)
    A[2:(D+1),2:(D+1)] <- C             
    
    reor <- c(1,varobs+1,dep+1)     
    A <- A[reor,reor]              
    A <- sweep.matrix(A,1:(q+1))             
    
    B <- A[1:(q+1),(q+2):(D+1)]       
    CR <- A[(q+2):(D+1),(q+2):(D+1)]       
    
    return(list(betas=B,resid=CR))
  }
  
  inv.alr <- function(x,pos){
    
    ad<-1/(rowSums(exp(x))+1)
    ax<-exp(x)*ad
    if(pos==1) {
      a<-cbind(ad,ax)
    }
    else { 
      if (dim(x)[2] < pos){
        a<-cbind(ax,ad)
      }   
      else {
        a<-cbind(ax[,1:(pos-1)],ad,ax[,pos:(dim(x)[2])])
      }
    }
    return(a)
  }
  
  ilr <- function (x){
    
    D <- length(x)
    z <- vector(mode="double",length=D-1)
    
    for (i in 1:(D-1)){
      z[i] <- sqrt((D-i)/(D-i+1))*log(x[i]/((prod(x[(i+1):D]))^(1/(D-i))))
    }
    return(as.numeric(z))
  }
  
  inv.ilr <- function(z){
    
    D <- length(z) + 1
    x <- vector(mode="double",length=D)
    
    x[1] <- exp(sqrt((D-1)/D)*z[1])
    for (i in 2:(D-1)){
      x[i] <- exp(-sum((1/sqrt((D-1:(i-1)+1)*(D-1:(i-1))))*z[1:(i-1)]) + sqrt(D-i)/sqrt(D-i+1)*z[i]) 
    }
    x[D] <- exp(-sum((1/sqrt((D-1:(D-1)+1)*(D-1:(D-1))))*z[1:(D-1)]))
    
    x <- x/sum(x)
    return(x)
  }
  
  
  ## Preliminaries ----  
  
  X <- as.data.frame(X)
  nn <- nrow(X); D <- ncol(X)
  
  X[X==label] <- NA
  X <- as.data.frame(apply(X,2,as.numeric))
  c <- apply(X,1,sum,na.rm=TRUE)
  
  if (nrow(dl)==1) dl <- matrix(rep(1,nn),ncol=1)%*%dl
  
  # Check for closure
  closed <- 0
  if (all( abs(c - mean(c)) < .Machine$double.eps^0.5 )) closed <- 1
  
  misspat <- as.data.frame(is.na(X)*1)
  misspat <- as.factor(do.call(paste,c(misspat,sep="")))
  levels(misspat) <- 1:(length(levels(misspat)))
  
  ## Ordinary lrEM  ----
  
  if (rob==FALSE){
    
    pos <- which(!is.na(colSums(X)))[1]
    if (is.na(pos)) stop("lrEM based on alr requires at least one complete column")
    
    cpoints <- log(dl)-log(X[,pos])-.Machine$double.eps
    cpoints <- cpoints[,-pos]
    
    X_alr <- log(X)-log(X[,pos]); X_alr <- as.matrix(X_alr[,-pos])
    tol <- tolerance
    nn <- nrow(X_alr); D <- ncol(X_alr)
    
    if (ini.cov != "multRepl"){
      if (inherits(try(solve(cov(X,use=ini.cov)),silent=TRUE),"try-error"))
        stop("ini.cov: too few complete cases for using 'complete.obs'")
      M <- matrix(colMeans(X_alr,na.rm=T),ncol=1)
      C <- cov(X_alr,use=ini.cov)}
    else {X.mr <- multRepl(X,label=NA,dl=dl,delta=delta)
          X.mr_alr <- t(apply(X.mr,1,function(x) log(x)-log(x[pos])))[,-pos]
          M <- matrix(colMeans(X.mr_alr,na.rm=T),ncol=1)
          C <- cov(X.mr_alr)}  
    
    iter_again <- 1
    niters <- 0
    alt.in <- FALSE
    alt.pat <- 0
    alt.mr <- 0
    
    while (iter_again == 1){
      
      niters <- niters + 1
      Mnew <- M                       
      Cnew <- C
      Y <- X_alr                              
      v <- matrix(0,D,D)
      
      for (npat in 2:length(levels(misspat))){                    
        i <- which(misspat==npat) 
        varobs <- which(!is.na(X_alr[i[1],]))
        varmiss <- which(is.na(X_alr[i[1],]))
        if (length(varobs) == 0){
          alt.in <- TRUE
          temp <- multRepl(X[i,,drop=FALSE],label=NA,dl=dl[i,,drop=FALSE],delta=delta)
          Y[i,] <- t(apply(temp,1,function(x) log(x)-log(x[pos])))[,-pos]
          if (niters == 1){
            alt.pat <- c(alt.pat,npat)
            alt.mr <- list(alt.mr,i)
          }
          break
        }
        sigmas <- matrix(0,ncol=D)
        B <- matrix(lm.sweep(M,C,varobs)[[1]],ncol=length(varmiss))
        CR <- lm.sweep(M,C,varobs)[[2]]
        Y[i,varmiss] <- matrix(1,nrow=length(i))%*%B[1,] + X_alr[i, varobs, drop=FALSE]%*%B[2:(length(varobs)+1),]
        sigmas[varmiss] <- sqrt(diag(as.matrix(CR)))
        for (j in 1:length(varmiss)){                                
          sigma <- sigmas[varmiss[j]]
          fdN01 <- dnorm((cpoints[i,varmiss[j]]-Y[i,varmiss[j]])/sigma)
          fdistN01 <- pnorm((cpoints[i,varmiss[j]]-Y[i,varmiss[j]])/sigma)
          Y[i,varmiss[j]] <- Y[i,varmiss[j]]-sigma*(fdN01/fdistN01)
        }
        v[varmiss,varmiss] <- v[varmiss,varmiss] + CR*length(i)
      }
      M <- matrix(colMeans(Y),ncol=1)
      dif <- Y - matrix(1,nrow=nn)%*%t(M)              
      PC <- t(dif)%*%dif                       
      C <- (PC+v)/(nn-1)   
      
      # Convergence check
      Mdif <- max(abs(M-Mnew))    
      Cdif <- max(max(abs(C-Cnew)))  
      if ((max(c(Mdif,Cdif)) < tol) | (niters == max.iter)) iter_again <- 0
    }
    
    Y <- inv.alr(Y,pos)
    
    for (i in 1:nn){
      if (any(is.na(X[i,]))){
        vbdl <- which(is.na(X[i,]))
        X[i,vbdl] <- (X[i,pos]/Y[i,pos])*Y[i,vbdl]
      }
    }     
  } # End ordinary lrEM
  
  ## Robust lrEM ----
  
  if (rob==TRUE){
    
    if (ini.cov == "multRepl") X.mr <- multRepl(X,label=NA,dl=dl,delta=delta)
    
    miss <- by(X,misspat,function(x) which(is.na(x[1,])))
    obs <- by(X,misspat,function(x) which(!is.na(x[1,])))
    
    iter_again <- 1
    niters <- 0
    X.old <- X
    alt.in <- FALSE
    alt.pat <- 0
    alt.mr <- 0
    
    nnn <- 0
    
    while (iter_again == 1){
      
      niters <- niters+1
      if (niters > 1) {X.old <- X; C.old <- C}
      
      for (npat in 2:length(levels(misspat))){                    
        if ((length(obs[[npat]]) == 1) & (!any(npat==alt.pat))){
          alt.in <- TRUE
          X[misspat==npat,] <- multRepl(X.old[misspat==npat,,drop=FALSE],label=NA,dl=dl[misspat==npat,,drop=FALSE],delta=delta)    
          alt.pat <- c(alt.pat,npat)
          alt.mr <- list(alt.mr,which(misspat==npat))
        }
        if (length(obs[[npat]]) > 1) {
          feeder <- X.old[,obs[[npat]]]
          for (m in 1:length(miss[[npat]])){
            p <- miss[[npat]][m]
            target <- X.old[,p]
            phi <- t(apply(cbind(dl=dl[misspat==npat,p],feeder[misspat==npat,]),1,ilr))
            regbasis <- as.data.frame(t(apply(cbind(target,feeder),1,ilr)))
            
            if (niters == 1){
              if (ini.cov == "complete.obs"){
                if (nrow(regbasis[misspat==1,]) > ncol(regbasis[misspat==1,]))
                  robreg <- rlm(V1 ~ .,data=regbasis[misspat==1,],method="MM",maxit = rlm.maxit)
                else
                  stop("ini.cov: too few complete cases. Use ini.cov = 'multRepl' instead")
              }
              if (ini.cov == "multRepl"){
                target <- X.mr[,p]
                feeder <- X.mr[,obs[[npat]]]
                regbasis.mr <- as.data.frame(t(apply(cbind(target,feeder),1,ilr)))
                robreg <- rlm(V1 ~ .,data=regbasis.mr,method="MM",maxit = rlm.maxit)
              }  
            }
            else
              robreg <- rlm(V1 ~ .,data=regbasis,method="MM",maxit = rlm.maxit)
            
            B <- matrix(robreg$coefficients,ncol=1)
            sigma <- robreg$s
            est <- cbind(V1=B[1,] + as.matrix(regbasis[misspat==npat,-1])%*%B[-1,],regbasis[misspat==npat,-1])
            est[,1] <- est[,1] - sigma*(dnorm((phi[,1]-est[,1])/sigma)/pnorm((phi[,1]-est[,1])/sigma))
            est <- t(apply(est,1,inv.ilr))
            est <- est[,1]*(feeder[misspat==npat,1]/est[,2])
            X[misspat==npat,p] <- est
          }
        }    
      }
      
      C <- cov(t(apply(X,1,ilr)))
      
      # Convergence check
      
      if (niters > 1)
        if((norm(C-C.old,type="F") < tolerance) | (niters == max.iter)) iter_again <- 0
      
    } 
  } # End robust lrEM
  
  ## Final section ----
  
  if (closed==1){
    X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
  }
  
  if (suppress.print==FALSE){
    if (alt.in) {
      cat("Warning: samples with only one observed component were found \n")
      for (i in 2:length(alt.pat)){
        cat(paste("  Pattern no.",alt.pat[i],"was imputed using multiplicative simple replacement \n"))
        cat("   Affected samples id: "); cat(alt.mr[[i]]); cat("\n\n")
      }
    }
  cat(paste("No. iterations to converge: ",niters,"\n\n"))
  }
  
  return(as.data.frame(X))  
  
}