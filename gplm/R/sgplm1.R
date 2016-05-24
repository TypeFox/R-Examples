"sgplm1" <- function(x,t,y,spar=NULL,df=4,family="gaussian",link="identity",
                     b.start=NULL,m.start=NULL,grid=NULL,
                     offset=0,method="speckman",weights=1,
                     weights.trim=1,weights.conv=1,max.iter=25,eps.conv=1e-8,
                     verbose=FALSE,...){
  
  x <- as.matrix(x)
  t <- as.matrix(t)
  y <- as.matrix(y)
  
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(t)
  if(q!=1) stop("sgplm1: t needs to be 1-dimensional")
  
  if (length(weights)==1){ weights <- rep(weights,n) }
  if (length(weights.trim)==1){ weights.trim <- rep(weights.trim,n) }
  if (length(weights.conv)==1){ weights.conv <- rep(weights.conv,n) }
  if (length(offset)==1){ offset<-rep(offset,n) }

  n.grid <- 0
  m.grid <- NULL
  havegrid<- !is.null(grid)
  if (havegrid) {
    grid <- as.matrix(grid)
    n.grid <- nrow(grid)
  }

  eps.conv <- eps.conv
  one <- rep(1,n)
  xnew <- matrix(rep(NA,n*p),nrow=n)

  if (!is.null(spar)){
    if (length(spar) < q){
      spar <- rep(spar,q)
    }
  }
  if (length(df) < q){
    df <- rep(df,q)
  }

  
  if (missing(b.start)) { b.start<- rep(0,p) }
  if (missing(m.start)) {
    if (family!="bernoulli"){
      m.start <- rep(glm.inverse.link(mean(y),family=family,link=link),n)
    }else{
      m.start <- rep(glm.inverse.link((2*y+1)/4),family=family,link=link)
    }
  }
  it <- 0
  stop.crit <- FALSE
  xb <- x %*% b.start + offset
  b <- b.start
  rss.start <- Inf

  while ( (stop.crit==FALSE)&(it< max.iter)){
    it <- it +1
    if (verbose){
      print( paste("sgplm1: iteration no.",as.character(it)) )
    }
    
    ll <- glm.lld(xb+m.start,y,family=family,link=link)
    zm <- m.start-ll$ll1.2 ## ll$ll1/ll$ll2
    z  <- xb+zm
    wnew <- as.vector(weights* ll$ll2)

    m.fit <- smooth.spline(x=t, y=zm, w=-wnew, spar=spar, df=df, ...)
    df.S <- m.fit$df
    m <- predict(m.fit,t)$y 

    for (j in 1:p){
      tmp <- smooth.spline(x=t, y=x[,j], w=-wnew, spar=spar, df=df, ...)
      xnew[,j] <- x[,j] - predict(tmp,t)$y
    }
    tmp <- smooth.spline(x=t, y=z, w=-wnew, spar=spar, df=df, ...)
    znew <- z - predict(tmp,t)$y
        
    if (method=="backfitting"){ 
      wnew<- x*as.vector(weights.trim*weights*ll$ll2)
    }
    else{
      wnew<- xnew*as.vector(weights.trim*weights*ll$ll2)
    }
    
    B  <- t(wnew) %*% xnew   ##; print(B)
    bv <- chol2inv(chol(-B))
    b  <- solve(B,t(wnew) %*% znew)

    eta <- xb+m
    mu  <- glm.link(eta, family=family, link=link)
    dev <- 2*sum(weights *glm.ll(y,y,family=family)-weights*glm.ll(mu,y,family=family))
    chgd <- abs((dev-rss.start)/rss.start)
    db   <- b-b.start
    chgb <- sqrt( sum(db*db)/sum(b.start*b.start) )
    dm   <- m-m.start
    chgm <- sqrt( sum(dm*dm)/sum(m.start*m.start) )

    if (it==1){
      stop.crit<- FALSE
    }
    else{
      stop.crit <- ( ((chgb<eps.conv)&(chgm<eps.conv)) | (chgd<eps.conv) )
    }
    if (verbose){
      print(paste("  Deviance: ",as.character(dev)))
      print(paste("  Change in b: ",as.character(chgb)))
      print(paste("  Change in m: ",as.character(chgm)))
      print(paste("  Change in deviance: ",as.character(chgd)))
    }
    b.start <- b
    m.start <- m
    rss.start <- dev
    xb<- x %*% b.start+offset
  }

  ## compute m on grid
  if (havegrid){
    m.grid <- predict(m.fit,grid)$y
  }

  ## compute df (residuals)
  xxnew <- xnew
  for (j in 1:p){
    tmp <- smooth.spline(x=t, y=xnew[,j],
                         w=-as.vector(weights* ll$ll2),spar=spar, df=df, ...)
    xxnew[,j] <- xnew[,j] - predict(tmp,t)$y
  }
  tmp  <- t(wnew) %*% (xnew-xxnew)
  df.res <- n + sum(diag(bv %*% tmp)) - df.S
  aic <- 2*(n-df)+dev

  return(list(b=b,bv=bv,m=m,m.grid=m.grid,it=it,df.residual=df.res,deviance=dev,aic=aic))
}

