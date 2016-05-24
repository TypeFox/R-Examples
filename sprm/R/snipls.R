snipls <-
function(data,eta,a,print=TRUE){
  Z <- data
  Xh <- scale(as.matrix(Z[,2:ncol(Z)]),center=TRUE,scale=FALSE)
  X0 <- Xh
  yh <- scale(as.vector(Z[,1]),center=TRUE,scale=FALSE)
  my <- attr(yh,"scaled:center")
  y0 <- yh
  Tpls <- NULL 
  W <- NULL 
  P <- NULL 
  C <- NULL 
  B <- NULL 
  bh <- 0
  Xev <- matrix(0,nrow=1,ncol=a)
  Yev <- matrix(0,nrow=1,ncol=a)
  oldgoodies <- NULL
  vars <- vector("list",2*a)
  for(i in 1:a){
    wh <- t(Xh)%*%yh
    wh <- wh/norm(wh,"F")
    goodies <- abs(wh)-eta*max(abs(wh))[1]
    wh <- goodies*sign(wh)
    goodies <- which((goodies>=0))
    goodies <- union(oldgoodies,goodies)
    oldgoodies <- goodies
    wh[setdiff(1:ncol(X0),goodies)] <- 0 
    th <- Xh%*%wh
    nth <- norm(th,"F")
    ch <- t(yh)%*%th/(nth^2)
    ph <- t(Xh)%*%th/(nth^2)
    ph[setdiff(1:ncol(X0),goodies)] <- 0 
    yh <- yh - th * as.numeric(ch)
	Xh <- Xh - th%*%t(ph)
    W <- cbind(W,wh)
    P <- cbind(P,ph) 
    C <- rbind(C,ch) 
    Tpls <- cbind(Tpls,th)
    Xev[i] <- (nth^2*norm(ph,"F")^2)/sum(X0^2)*100
    Yev[i] <- sum(nth^2*as.numeric(ch^2))/sum(y0^2)*100
    if(print==TRUE){cat("Variables retained for ",i," latent variable(s):","\n",
                     colnames(X0)[goodies],".\n")}
    vars[[2*(i-1)+1]] <- colnames(X0)[goodies] # paste(colnames(X0)[goodies],collapse=", ")
    vars[[2*i]] <- goodies
  }
  if(length(goodies)>0){
    R <- W %*% solve(t(P)%*%W)
    B <- R%*%C
  } else{
    B <- matrix(0,nrow=ncol(Xh),ncol=1)
    R <- B
    Tpls <- matrix(0,nrow=nrow(Xh),ncol=a)
  }
  yp <- X0%*%B + my
  if (any(is.nan(Tpls))){
	  stop("NaN generated in Tpls")
  }
  if(length(vars[[2*a]])==0){
    stop("No variables have been retained in Sparse PRM model!")
  }
  return(list(W=W,loadings=P, C=C, scores=Tpls, coefficients=B, Xev=Xev, Yev=Yev, Vars=vars, fitted.values=yp, R=R ))
}
