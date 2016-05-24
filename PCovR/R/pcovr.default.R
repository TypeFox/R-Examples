pcovr.default <-
  function(X,Y,modsel="seq",Rmin=1,Rmax=ncol(X)/3,R=NULL,weight=NULL,rot="varimax", target=NULL, prepX="stand",prepY="stand", ratio=ErrorRatio(X,Y,Rmin,Rmax,prepX,prepY), fold="LeaveOneOut",zeroloads=ncol(X)){
    J <- ncol(X)
    N <- nrow(X)
    K <- ncol(Y)
    if (is.null(ncol(Y))==T){
      K <- 1
    }
    Jlabel <- colnames(X)
    Klabel <- colnames(Y)
    
    if (is.null(R)==F){
      Rmin <- Rmax <- R
    }
    vec <- Rmin:Rmax
    
    
    # PREPROCESSING
    X <- switch(prepX,
                stand=nrm2(scale(X, center=T, scale=F))*N^(1/2),
                cent=scale(X, center=T, scale=F))
    Y <- switch(prepY,
                stand=nrm2(scale(Y, center=T, scale=F))*N^(1/2),
                cent=scale(Y, center=T, scale=F))
    Y <- array(Y,c(N,K))
    
    
    # MODEL SELECTION
    AlphaMaxLik <- (SUM(X)$ssq/ (SUM(X)$ssq + SUM(Y)$ssq*ratio))
    if (length(weight)==1){
      a <- weight
    } else {
      a <- AlphaMaxLik
    }
    
    if (modsel=="sim"){
      if (is.null(weight)){
        a <- seq(from=.01, to=1, by=.01)
      } else {
        a <- weight
      }
      CV <- array(NA,c(length(a),length(vec)))
      for (l in 1:length(vec)){
        for (w in 1:length(a)){
          para <- pcovr_est(X,Y, r = vec[l], a = a[w], cross = TRUE,fold)
          CV[w,l] <- para$Qy2
        }
      }
      alpha <- a[arrayInd(which.max(CV),dim(CV))[1]]
      R <- vec[arrayInd(which.max(CV),dim(CV))[2]]
      
    } else if (modsel=="seqRcv"){
      if (length(weight)>1){
        print("WARNING: Weighting values overruled by maximum likelihood procedure")
      }
      CV <- array(NA,c(1,length(vec)))
      for (l in 1:length(vec)){
        para <- pcovr_est(X,Y,vec[l],a,TRUE,fold)
        CV[1,l] <- para$Qy2
      }
      R <- vec[which.max(CV)]
      alpha <- a
      
    } else if (modsel=="seq") {
      if (length(weight)>1){
        print("WARNING: Weighting values overruled by maximum likelihood procedure")
      }
      vec <- c(vec[1]-1,vec,vec[length(vec)]+1)
      VAF <- matrix(0,1,length(vec))
      scr <- array(NA,c(1,length(vec)))
      for (l in 1:length(vec)){
        if (vec[l]>0){
          para <- pcovr_est(X,Y,vec[l], a = a, cross = FALSE,fold)
          VAF[1,l] <- a*para$Rx2 + (1-a)*para$Ry2
        }
      }
      for (u in 2:(length(vec)-1)){
        scr[,u]=(VAF[u]-VAF[u-1])/(VAF[u+1]-VAF[u])
      }
      R <- vec[which.max(scr)]
      alpha <- a
    } else if (modsel=="seqAcv"){
      if (is.null(weight)){
        a <- seq(from=.01, to=1, by=.01)
      } else {
        a <- weight
      }
      vec <- c(vec[1]-1,vec,vec[length(vec)]+1)
      VAF <- matrix(0,1,length(vec))
      scr <- array(NA,c(1,length(vec)))
      CV <- array(NA,c(length(a),1))
      for (l in 1:length(vec)){
        if (vec[l]>0){
          para <- pcovr_est(X,Y,vec[l], a = AlphaMaxLik, cross = FALSE,fold)
          VAF[1,l] <- AlphaMaxLik*para$Rx2 + (1-AlphaMaxLik)*para$Ry2
        }
      }
      for (u in 2:(length(vec)-1)){
        scr[,u]=(VAF[u]-VAF[u-1])/(VAF[u+1]-VAF[u])
      }
      R <- vec[which.max(scr)]
      for (w in 1:length(a)){
        para <- pcovr_est(X,Y,R,a[w],TRUE,fold)
        CV[w,1] <- para$Qy2
      }
      alpha <- a[which.max(CV)]
    } else {
      print("ERROR: This model selection procedure does not exist")
    }
    
    
    # ANALYSES
    results <- list()
    para <- pcovr_est(X,Y,r = R, a = alpha, cross = FALSE,fold)
    Te <- para$Te * N^(1/2)
    W <- para$W * N^(1/2)
    Px <- para$Px * N^(-1/2)
    Py <- para$Py * N^(-1/2)
    
    #rotate parameters
    if ((R>1) & (rot!="none")){
      Bpx <- t(Px)
      if (rot=="quartimin"){
        rotation <- GPFoblq(Bpx, Tmat=diag(ncol(Bpx)))
      }
      if (rot=="varimax"){
        rotation <- GPForth(Bpx, Tmat=diag(ncol(Bpx)))
      }
      if (rot=="targetT"){
        rotation <- targetT(Bpx, Tmat=diag(ncol(Bpx)), t(target))
      }
      if (rot=="targetQ"){
        rotation <- targetQ(Bpx, Tmat=diag(ncol(Bpx)), t(target))
      }
      if (rot=="simplimax"){
        rotation <- simplimax(Bpx, Tmat=diag(ncol(Bpx)),k=zeroloads)
      }
      if (rot=="wvarim"){
        rotation <- wvarim(Bpx)
      }
      if (rot=="promin"){
        rotation <- promin(Bpx)
      }
      rotE <- t(ginv(rotation$Th))
      Px <- t(rotation$loadings)
      Te <- Te %*% t(ginv(rotE))
      W <- W %*% t(ginv(rotE))
      Py <- t(rotE) %*% Py
    }
    
    #results
    Rx2 <- para$Rx2
    Ry2 <- para$Ry2
    Rstr <- paste("component",1:R,sep="") 
    Px <- data.frame(Px)
    rownames(Px) <- Rstr
    colnames(Px) <- Jlabel
    Py <- data.frame(Py)
    rownames(Py) <- Rstr
    colnames(Py) <- Klabel
    Te <- data.frame(Te)
    colnames(Te) <- Rstr
    W <- data.frame(W)
    colnames(W) <- Rstr
    rownames(W) <- Jlabel
    
    if (modsel == "sim"){
      CV <- data.frame(CV)
      rownames(CV) <- paste("alpha=",a,": ",sep="")
      colnames(CV) <- paste(vec,"component(s)",sep="")
      results <- list(Px=Px,Py=Py,Te=Te,W=W,Rx2=Rx2,Ry2=Ry2,Qy2=CV,alpha=alpha,R=R,modsel=modsel,rot=rot,prepX=prepX,prepY=prepY,Rvalues=vec,Alphavalues=a)
    } else if (modsel=="seqRcv"){
      CV <- data.frame(CV)
      rownames(CV) <- paste("alpha=",alpha,": ",sep="")
      colnames(CV) <- paste(vec,"component(s)",sep="")
      results <- list(Px=Px,Py=Py,Te=Te,W=W,Rx2=Rx2,Ry2=Ry2,Qy2=CV,alpha=alpha,R=R,modsel=modsel,rot=rot,prepX=prepX,prepY=prepY,Rvalues=vec,Alphavalues=a)
    } else if (modsel=="seqAcv"){
      CV <- data.frame(CV)
      rownames(CV) <- paste("alpha=",a,": ",sep="")
      colnames(CV) <- paste(R,"component(s)",sep="")
      VAF <- data.frame(VAF)
      rownames(VAF) <- paste("alpha=",a[which.min(abs(a - AlphaMaxLik))],": ",sep="")
      colnames(VAF) <- paste(vec,"component(s)",sep="")
      results <- list(Px=Px,Py=Py,Te=Te,W=W,Rx2=Rx2,Ry2=Ry2,Qy2=CV,VAFsum=VAF,alpha=alpha,R=R,modsel=modsel,rot=rot,prepX=prepX,prepY=prepY,Rvalues=vec,Alphavalues=a)
    } else {
      VAF <- data.frame(VAF)
      rownames(VAF) <- paste("alpha=",alpha,": ",sep="")
      colnames(VAF) <- paste(vec,"component(s)",sep="")
      results <- list(Px=Px,Py=Py,Te=Te,W=W,Rx2=Rx2,Ry2=Ry2,VAFsum=VAF,alpha=alpha,R=R,modsel=modsel,rot=rot,prepX=prepX,prepY=prepY,Rvalues=vec,Alphavalues=a)
    }
    class(results) <- "pcovr"
    return(results)
  }