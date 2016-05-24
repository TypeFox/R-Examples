lrDA <-
  function(X,label=NULL,dl=NULL,ini.cov=c("lrEM","complete.obs","multRepl"),delta=0.65,n.iters=1000,m=1,store.mi=FALSE){
    
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
    
    if ((store.mi==TRUE) & (m==1)) store.mi <- FALSE
    
    lm.sweep <- function(M,C,varobs){
      
      sweep.matrix <- function(A,ind){ 
        
        nn <- nrow(A); p <- ncol(A)
        S <- A
        
        for (j in ind){                       
          S[j,j] <- -1/A[j,j]     
          for (i in 1:p) {           
            if (i != j){
              S[i,j] <- -A[i,j]*S[j,j]
              S[j,i] <- S[i,j]
            }
          }
          for (i in 1:p){           
            if (i != j){
              for (k  in 1:p){
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
      
      p <- length(M)                      
      q <- length(varobs)                     
      i <- rep(1,p)                  
      i[varobs] <- i[varobs]-1
      dep <- which(i!=0)
      ndep <- length(dep)                 
      
      A <- matrix(0,p+1,p+1)             
      A[1,1] <- -1                    
      A[1,2:(p+1)] <- M                
      A[2:(p+1),1] <- matrix(M,ncol=1)
      A[2:(p+1),2:(p+1)] <- C             
      
      reor <- c(1,varobs+1,dep+1)     
      A <- A[reor,reor]              
      A <- sweep.matrix(A,1:(q+1))             
      
      B <- A[1:(q+1),(q+2):(p+1)]       
      CR <- A[(q+2):(p+1),(q+2):(p+1)]       
      
      return(list(betas=B,resid=CR))
    }
    
    inv.raw <- function(Y,X,pos,closed,nn,c){
      
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
      
      Y <- inv.alr(Y,pos)
      
      for (i in 1:nn){
        if (any(is.na(X[i,]))){
          vbdl <- which(is.na(X[i,]))
          X[i,vbdl] <- (X[i,pos]/Y[i,pos])*Y[i,vbdl]
        }
      }    
      if (closed==1){
        X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
      }
      return(as.data.frame(X))
    }
    
    ini.cov <- match.arg(ini.cov)
    
    X <- as.data.frame(X)
    nn <- nrow(X); p <- ncol(X)
    
    X[X==label] <- NA
    X <- apply(X,2,as.numeric)
    c <- apply(X,1,sum,na.rm=TRUE)
    
    if (nrow(dl)==1) dl <- matrix(rep(1,nn),ncol=1)%*%dl
    
    # Check for closure
    closed <- 0
    if (all( abs(c - mean(c)) < .Machine$double.eps^0.5 )) closed <- 1
    
    pos <- which(!is.na(colSums(X)))[1]
    if (is.na(pos)) stop("lrDA requires at least one fully observed column")
    
    cpoints <- log(dl)-log(X[,pos])-.Machine$double.eps
    cpoints <- cpoints[,-pos]
    
    X_alr <- log(X)-log(X[,pos]); X_alr <- as.matrix(X_alr[,-pos])
    nn <- nrow(X_alr); p <- ncol(X_alr)
    
    if (ini.cov == "complete.obs"){
      if (inherits(try(solve(cov(X,use=ini.cov)),silent=TRUE),"try-error"))
        stop("ini.cov: too few complete cases for using 'complete.obs'")
      M <- matrix(colMeans(X_alr,na.rm=T),ncol=1)
      C <- cov(X_alr,use=ini.cov)}
    if (ini.cov == "multRepl"){
      X.mr <- multRepl(X,label=NA,dl=dl,delta=delta)
      X.mr_alr <- t(apply(X.mr,1,function(x) log(x)-log(x[pos])))[,-pos]
      M <- matrix(colMeans(X.mr_alr,na.rm=T),ncol=1)
      C <- cov(X.mr_alr)}
    if (ini.cov == "lrEM"){
      X.em <- lrEM(X,label=NA,dl=dl,ini.cov="multRepl",delta=delta,suppress.print=TRUE)
      X.em_alr <- t(apply(X.em,1,function(x) log(x)-log(x[pos])))[,-pos]
      M <- matrix(colMeans(X.em_alr,na.rm=T),ncol=1)
      C <- cov(X.em_alr)}
    
    misspat <- as.data.frame(is.na(X)*1)
    misspat <- as.factor(do.call(paste,c(misspat,sep="")))
    levels(misspat) <- 1:(length(levels(misspat)))
    
    t <- 0
    k <- 0
    runs <- 0
    alt.in <- FALSE
    alt.pat <- 0
    alt.mr <- 0
    
    if (m > 1){
      imputed <- matrix(0,nrow=m,ncol=sum(is.na(X_alr)))
      if (store.mi==TRUE) mi.list <- vector(mode="list",m)
    }
    
    while (t <= n.iters*m){
      
      Y <- X_alr                              
      runs <- runs + 1
      
      # I-step
      
      for (npat in 2:length(levels(misspat))){                     
        i <- which(misspat==npat) 
        varobs <- which(!is.na(X_alr[i[1],]))
        varmiss <- which(is.na(X_alr[i[1],]))
        if (length(varobs) == 0){
          alt.in <- TRUE
          temp <- multRepl(X[i,,drop=FALSE],label=NA,dl=dl[i,,drop=FALSE],delta=delta)
          Y[i,] <- t(apply(temp,1,function(x) log(x)-log(x[pos])))[,-pos]
          if (runs == 1){
            alt.pat <- c(alt.pat,npat)
            alt.mr <- list(alt.mr,i)
          }
          break
        }
        sigmas <- matrix(0,ncol=p)
        B <- matrix(lm.sweep(M,C,varobs)[[1]],ncol=length(varmiss))
        CR <- lm.sweep(M,C,varobs)[[2]]
        Y[i,varmiss] <- matrix(1,nrow=length(i))%*%B[1,] + X_alr[i, varobs, drop=FALSE]%*%B[2:(length(varobs)+1),]
        sigmas[varmiss] <- sqrt(diag(as.matrix(CR)))
        for (j in 1:length(varmiss)){                                
          sigma <- sigmas[varmiss[j]]
          Y[i,varmiss[j]] <- rtruncnorm(1,-Inf,cpoints[i,varmiss[j]],Y[i,varmiss[j]],sigma)
        }
      }
      
      if ((t%in%((1:m)*n.iters)) & (m > 1)){
        k <- k + 1
        imputed[k,] <- Y[which(is.na(X_alr))]
        if (store.mi==TRUE){
         mi.list[[k]] <- Y
        }
      }

      # P-step
      
      C <- rWishart(1,nn-1,solve(nn*cov(Y)))[,,1]
      M <- mvrnorm(1,colMeans(Y),1/nn*C)

    t <- t + 1

    }
    
    if ((m > 1) & (store.mi == FALSE)) Y[which(is.na(X_alr))] <- colMeans(imputed) # MI estimates
    
    if (store.mi==FALSE) X <- inv.raw(Y,X,pos,closed,nn,c)
    
    if (store.mi==TRUE) X <- lapply(mi.list,FUN=function(x) inv.raw(x,X,pos,closed,nn,c))
    
    if (alt.in) {
      cat("Warning: samples with only one observed component were found \n")
      for (i in 2:length(alt.pat)){
        cat(paste("  Pattern no.",alt.pat[i],"was imputed using multiplicative simple replacement \n"))
        cat("   Affected samples id: "); cat(alt.mr[[i]]); cat("\n\n")
      }
    }
    return(X)
  }
