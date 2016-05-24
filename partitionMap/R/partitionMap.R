partitionMap <-
function(X,Y,XTEST=NULL,YTEST=NULL,method="pm",dimen=2,force=TRUE,ntree=100,plottrain=TRUE,addjitter=0.03,...){

  if(class(Y)!="factor"){
    warning("\n converting response to a factor variable")
    Y <- as.factor(Y)
  }
    
  n <- nrow(X)
  if(!is.null(XTEST)) nall <- nrow(X)+nrow(XTEST) else nall <- n
  p <- ncol(X)


  rfall <- randomForest(X, as.factor(as.character(Y)) ,ntree=100,keep.forest=TRUE,keep.inbag=TRUE,...)

  predall <- predict(rfall,newdat=if(!is.null(XTEST)) rbind(X,XTEST) else X,nodes=TRUE,predict.all=TRUE)
  NO <- attr(predall,"nodes")
  Z <- matrix(0,nrow=nall,ncol= sum(apply(NO,2, function(x) length(unique(x)))))
  nc <- 0
  for (treec in 1:rfall$ntree){
    for (nod in unique(NO[,treec])){
      nc <- nc+1
      ind <- as.numeric(which(NO[,treec]==nod))
      Z[ind,nc] <- 1
    }
  }
  Z <- cbind(rep(1,nrow(Z)),Z)
  

  if(!is.null(XTEST)) ZTEST <- Z[-(1:n),]
  Z <- Z[1:n,]
  q <- ncol(Z)
  ZN <- sweep(Z,1,apply(Z,1,sum),FUN="/")
  

  dias <- apply(Z,1,sum)
  diar <- apply(Z,2,sum)

  A2 <- sweep( sweep( Z,1,sqrt(dias),FUN="/") ,  2, sqrt(diar),FUN="/")
  cat("\n\n ******")

  SAMPLES <- scale(diag(1/sqrt(dias)) %*% svd(A2 + matrix(0.01*rnorm(nrow(A2)*ncol(A2)),nrow=nrow(A2)),nu=3,nv=0)$u[,1+(1:dimen)])
  cat("\n\n ******")

  SAMPLES <- SAMPLES[,1:dimen]
  RULES <- sweep(  t(Z), 1, apply(t(Z),1,sum), FUN="/") %*% SAMPLES
  SAMPLES <- sweep( Z,1,apply(Z,1,sum),FUN="/") %*% RULES
  if(!is.null(XTEST)) SAMPLESTEST <- sweep( ZTEST,1,apply(ZTEST,1,sum),FUN="/") %*% RULES

      cat("\n\n ******")

  sdvec <- apply(SAMPLES,2,sd)
  SAMPLES <- (SAMPLES)/max(sdvec)
  if(!is.null(XTEST)) SAMPLESTEST <- (SAMPLESTEST)/max(sdvec)
  RULES <- (RULES)/max(sdvec)
  
      cat("\n\n ******")

  if(method=="pm"){
    
    nsimloop <- 1000
    eps <- 0.01
    DIST <- (outer(Y,Y,FUN="!="))
    change <- 1
    simloop <- 0
    
    nY <- length(unique(Y))
    ZS <- matrix(0,nrow=nY,ncol=ncol(Z))
    for (i in 1:length(unique(Y))){
      ZS[i,] <- apply(Z[ Y==unique(Y)[i], ,drop=FALSE],2,sum)
    }
    ds <- diag(ZS %*% t(ZS))
    ZS <- sweep(ZS,1,apply(ZS,1,sum),FUN="/")

    cat("\n\n ******")
    cat("\n  start iterations... ")

    if(force){
      SAMPLESY0 <- matrix(0,nrow=nY,ncol=dimen)
      for (ii in 1:length(unique(Y))){
        SAMPLESY0[ii,] <- apply( SAMPLES[ unique(Y)[ii]==Y, ,drop=FALSE],2,mean) * 2
      }
      eps <- 0.01
      change <- 1
      simloop <- 0
      SAMPLESY <- SAMPLESY0
      DIST <- 1-diag(nY)
      
      while(change > 10^(-6)){
        simloop <- simloop+1
        B <- SAMPLESY %*% t(SAMPLESY)
        DCURR <- outer(diag(B),rep(1,nY)) + outer(rep(1,nY),diag(B)) -2* B
        grad <- matrix(nrow=nY,ncol=dimen)
        grade <- numeric(dimen)
        
        for (i in 1:nY){
          for (dimc in 1:dimen){
            grade[dimc] <- sum(DIST[i,]/pmax(eps,DCURR[i,])^2*2* (SAMPLESY[i,dimc] - SAMPLESY[,dimc]))
          }
          
          grad0 <- (SAMPLESY0[i,1:dimen]-SAMPLESY[i,1:dimen])
          
          grad[i,] <- grade + grad0
          
        }
        l2 <- sqrt(sum(apply(grad^2,1,sum)))
        SAMPLESY <- SAMPLESY + 1/(1+simloop/(nsimloop)) *0.1* grad/(0.00001+l2)
        
        B <- SAMPLESY %*% t(SAMPLESY)
        DCURR <- outer(diag(B),rep(1,nY)) + outer(rep(1,nY),diag(B)) -2* B
        
        if(round(simloop/10)==simloop/10 ){
          RULES <- t( sweep(ZS,2,apply(ZS,2,sum),FUN="/") ) %*% SAMPLESY
          
          change <- mean( (SAMPLESY0 - ZS %*% RULES )^2)/ mean( (SAMPLESY0 )^2)
          SAMPLESY0 <- ZS %*% RULES 
        }
        
      }
      RULES <- t( sweep(ZS,2,apply(ZS,2,sum),FUN="/") ) %*% SAMPLESY
      
    }else{
      
      B <-  ZS %*% sweep( t(ZS), 1,apply(ZS,2,sum),FUN="/")
      CENTER <- diag(nY) - outer(rep(1,nY),rep(1,nY),FUN="*")/nY
      SAMPLESY <-  eigen(CENTER %*% B %*% CENTER)$vec[,1:max(2,min(nY-1,3)),drop=FALSE]
      RULES <- t( sweep(ZS,2,apply(ZS,2,sum),FUN="/") ) %*% SAMPLESY
    }
    
    cat("\n  end iterations... ")
    SAMPLESY <- SAMPLESY[,1:dimen]
    RULES <- RULES[,1:dimen]
    
    SAMPLES <- sweep( Z,1,apply(Z,1,sum),FUN="/") %*% RULES
    if(!is.null(XTEST)) SAMPLESTEST <- sweep( ZTEST,1,apply(ZTEST,1,sum),FUN="/") %*% RULES

    sdvec <- apply(SAMPLES,2,sd)
    SAMPLESY <- (SAMPLESY)/max(sdvec)
    SAMPLES <- (SAMPLES)/max(sdvec)
    if(!is.null(XTEST))SAMPLESTEST <- (SAMPLESTEST)/max(sdvec)
    RULES <- (RULES)/max(sdvec)
  }
  ret <- if(!is.null(XTEST)) list(Samples=SAMPLES,Rules=RULES,Z=Z,Samplestest=SAMPLESTEST,Ztest=ZTEST)  else list(Samples=SAMPLES,Rules=RULES,Z=Z, rf =rfall)

  
  if(plottrain){
    par(mfrow=c(1,2 +as.numeric(!is.null(XTEST))))
    plot(ret$Samples + addjitter*sd(as.vector(ret$Samples))*matrix(rnorm(dimen*length(Y)),ncol=dimen)   ,col=Y,pch=20,cex=1.5,main="Training Data",
         xlab="Dimension 1",ylab="Dimension 2")
    points(ret$Rules,pch=".")
    if(!is.null(XTEST)){
      plot(ret$Samplestest +  addjitter*sd(as.vector(ret$Samplestest))*matrix(rnorm(dimen*nrow(XTEST)),ncol=dimen) ,col=if(!is.null(YTEST)) YTEST else "darkgrey",pch=20,cex=1.5,main="Test Data", xlab="Dimension 1",ylab="Dimension 2")
      points(ret$Rules,pch=".")
    }
    plot(ret$Samples,col=Y,pch=20,cex=1.5,xlab="",ylab="",type="n",axes=FALSE)
    legend(quantile(ret$Samples[,1],0),quantile(ret$Samples[,2],1),unique(Y),
           col=unique(Y),fill=unique(Y),border=0)
    par(mfrow=c(1,1))
  }
  
  return(ret)
}

