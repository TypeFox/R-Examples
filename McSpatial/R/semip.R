semip <- function(form,nonpar,conpar=NULL,window1=.25,window2=.25,bandwidth1=0,bandwidth2=0,kern="tcub",distance="Mahal",targetfull=NULL,
  print.summary=TRUE, data=NULL) {

  xmat <- model.frame(form,data=data)
  y <- xmat[,1]
  xmat <- as.matrix(xmat[,-1])
  nk = ncol(xmat)
  n = length(y)
  xname <- colnames(xmat)
  emat <- xmat
  basedata <- data.frame(y,xmat,model.frame(nonpar,data=data))
  if (!identical(conpar,NULL)) {basedata <- data.frame(basedata, model.frame(conpar,data=data)) }

  tvect <- targetfull
  if (!identical(targetfull,NULL)&!identical(targetfull,"alldata")){
    if (!identical(conpar,NULL)){tvect <- targetfull$obs}
    if (identical(conpar,NULL)){tvect <- targetfull$target}
  }


  if (identical(conpar,NULL)) {
    formy <- update(nonpar, y~., env=basedata)
    fit <- lwr(formy,window=window1,bandwidth=bandwidth1,kern=kern,distance=distance,target=tvect,data=basedata)
    ey <- y-fit$yhat
    for (j in seq(1:nk)) {
      basedata$x <- xmat[,j]
      formx <- update(nonpar, x~., env=basedata)
      fit <- lwr(formx,window=window1,bandwidth=bandwidth1,kern=kern,distance=distance,target=tvect,data=basedata)
      emat[,j] <- basedata$x-fit$yhat
    }
  }
  if (!identical(conpar,NULL)) {
    formy <- update(conpar, y~., env=basedata)
    fit <- cparlwr(formy,nonpar,window=window1,bandwidth=bandwidth1,kern=kern,distance=distance,targetobs=tvect,data=basedata)
    ey <- y-fit$yhat
    for (j in seq(1:nk)) {
      basedata$x <- xmat[,j]
      formx <- update(conpar, x~., env=basedata)
      fit <- cparlwr(formx,nonpar,window=window1,bandwidth=bandwidth1,kern=kern,distance=distance,targetobs=tvect,data=basedata)
      emat[,j] <- basedata$x-fit$yhat
    }
  }

  xbfit <- lm(ey~as.matrix(emat)+0)
  names(xbfit$coefficients) <- xname
  xcoef <- xbfit$coef
  xbhat <- as.vector(xmat%*%xbfit$coef)
  xx <- solve(crossprod(as.matrix(emat)))
  basedata$e <- y-xbhat

  if (identical(conpar,NULL)) {
    forme <- update(nonpar, e~., data=basedata)
    npfit <- lwr(forme,window=window1,bandwidth=bandwidth1,kern=kern,distance=distance,target=tvect,data=basedata)
  }
  if (!identical(conpar,NULL)) {
    forme <- update(conpar, e~., data=basedata)
    npfit <- cparlwr(forme,nonpar,window=window1,bandwidth=bandwidth1,kern=kern,distance=distance,targetobs=tvect,data=basedata)
  }
  nphat <- npfit$yhat
  
  df1 = npfit$df1 + nk
  df2 = npfit$df2 + nk

  rss = sum((y-xbhat-nphat)^2)
  sig2 = rss/(n-2*df1 + df2)
  gcv = n*rss/((n-df1)^2)
  vmat <- sig2*xx
  if (print.summary==TRUE) {
    semat <- sqrt(diag(vmat))
    outmat <- cbind(xcoef, semat, xcoef/semat, 2*(1-pnorm(abs(xcoef)/semat)) )
    rownames(outmat) <- names(xcoef)
    colnames(outmat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
    cat("Parametric Portion","\n")
    cat(" ", "\n")
    print(outmat)
  }

  out <- list(xcoef,vmat,xbhat,nphat,npfit$yhat.se,npfit,df1,df2,sig2,gcv)
  names(out) <- c("xcoef","vmat","xbhat","nphat","nphat.se","npfit","df1","df2","sig2","gcv")
  return(out)
}
 
