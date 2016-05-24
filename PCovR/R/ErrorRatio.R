ErrorRatio <-
  function(X,Y,Rmin=1,Rmax=ncol(X)/3,prepX="stand",prepY="stand"){
    J <- ncol(X)
    N <- nrow(X)
    
    #Preprocessing
    X <- switch(prepX,
                stand=nrm2(scale(X, center=T, scale=F))*N^(1/2),
                cent=scale(X, center=T, scale=F))
    Y <- switch(prepY,
                stand=nrm2(scale(Y, center=T, scale=F))*N^(1/2),
                cent=scale(Y, center=T, scale=F))
    K <- ncol(Y)
    Y <- array(Y,c(N,K))
    
    #Explained variance Y
    xnam <- paste("V", 1:J, sep="")
    if (K==1){
      ynam <- paste("V", (J+1), sep="")
      formula1 <- as.formula(paste(paste(ynam," ~ "), paste(xnam, collapse= "+")))    #Explained variance Y
    } else {
      ynam <- paste("V", (J+1):(J+K), sep="")
      formula1 <- as.formula(paste(paste("cbind(",paste(ynam, collapse =','),") ~ "), paste(xnam, collapse= "+")))    #Explained variance Y
    }
    data <- as.data.frame(cbind(X,Y))
    colnames(data) <- paste("V", 1:(J+K), sep="")
    reg <- lm(formula1, data)
    if (K==1){
      Ry2 <- SUM(as.matrix(X) %*% as.matrix(reg$coefficients[2:(J+1)]))$ssq / SUM(Y)$ssq
    } else {
      Ry2 <- SUM(as.matrix(X) %*% as.matrix(reg$coefficients[2:(J+1),]))$ssq / SUM(Y)$ssq
    }
    ery <- 1-Ry2
    
    #Explained variance X
    sing <- svd(X,nv=9)
    vec <- Rmin:Rmax
    vec <- c(vec[1]-1,vec,vec[length(vec)]+1)
    VAF <- c(0,cumsum(sing$d^2)/sum(sing$d^2))
    VAF <- VAF[vec+1]
    scr <- array(NA,c(1,length(vec)))
    for (u in 2:(length(vec)-1)){
      scr[,u]=(VAF[u]-VAF[u-1])/(VAF[u+1]-VAF[u])
    }
    R <- vec[which.max(scr)]
    erx <- 1-VAF[which.max(scr)]
    
    #error variance ratio
    ratio <- erx/ery
    return(ratio)
  }
