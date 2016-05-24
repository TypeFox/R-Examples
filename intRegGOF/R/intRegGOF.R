# compIntRegGOF.R
#   Interface R to computed integrated regression.
# NOTAS:
# ERRORES:






getModelFrame <- function(obj)
##
# obj = an object (lm, glm, nlm) with a model.frame.
{
  return( as.data.frame( model.frame(obj) ) )
}


getModelCovars <- function(obj)
##
# obj = an object (lm, glm, nlm) with a model.frame.
#   RETURN:
#     Names (strings) of all covariables used in model.
{
  # TODO: Para la clase nls no se borra la respuesta
  ft <- terms(getModelFrame(obj))
  ##  return( all.vars(delete.response(ft)) )
  return( all.vars(ft)[-1] )
}


getContVar <- function(df,vars=NULL)
##
# df = data frame or model fit object (lm, glm,nlm)
# vars = vector with variable names in df 
#   RETURN:
#     Subset of vars that are in data frame df that are continuous(numeric).
#   .
{
  if( is.data.frame(df) && is.null(vars) )
    stop("Error: getContVar(data.frame,NULL).")
  if(is.null(vars))
    vars <- getModelCovars(df)
  if(!is.data.frame(df))
    df <- eval( df$call$data,envir=parent.frame() )
  stopifnot(all(vars %in% names(df)))
  if (length(vars)==1) {
    if(is.numeric(df[,vars])) return(vars) else stop("No numeric covariables")
  } else
    return( vars[ sapply(df[,vars],is.numeric)==TRUE ] )
}


getResiduals <- function(obj,type="response")
##
# obj = an object (lm, glm, nlm) with a model.frame.
# type = type of residual, see 'glm' or 'lm'
#   RETURN:
#      Vector of residuals for model in obj.
{
  if( "glm" %in% class(obj) )
    hatEps <- residuals(obj,type="response")
  else
    hatEps <- residuals(obj)
  return( hatEps )
}


getModelWeights <- function(obj)
##
# obj = an object (lm, glm, nlm) with a model.frame.
#   RETURN:
#      Vector of weights for model in obj.
{
  mf <- getModelFrame(obj)
  weig <- model.weights(mf)
  if( is.null(weig) ) 
    weig <- rep(1,nrow(mf))
  return( weig )
}


compIntRegProc <- function(y,xord,weig=rep(1,length(y)))
##
# y = vector with y values.
# xord = a vector with the summation order or a list provided by getLessThan()
# weig = weigths associated with every y point.
{
  n <- length(y)
  stopifnot(n==length(xord),n==length(weig))
  if( is.list(xord) ) {
    res <- mvCumSum(y*weig,xord)/sqrt(n)
  } else{
    res <- cumsum(y[xord]*weig[xord])/sqrt(n)
  }
  return( res )
}


compBootSamp <- function(obj,datLT,B=999,LINMOD=FALSE)
# obj = an lm, glm or nlm fit object.
# datLT = a vector with the summation order or a list provided by getLessThan()
# B = bootstrap sample size
# LINMOD = when true Linear Model matrix fitting equation is used.
# return a matrix with colum names "Kn" and "Wn"
##  OBS: nlm is not checked
{
  hatEps <- getResiduals(obj)
  hatY <- fitted(obj)
  n <- length(hatY)
  weig <- getModelWeights(obj)
  bootSamp <- NULL
  if( LINMOD )
  {
    xMat <- model.matrix(obj)
    wMat <- diag(weig)
    hMat <- xMat %*% solve(t(xMat) %*% wMat %*% xMat) %*% t(xMat) %*% wMat
    for(i in 1:B)
    {
      astY <- hatY + hatEps * rWildBoot(n)
      hatAstY <- hMat %*% astY
      ##TODO: compIntRegProc(as.vector(astY-hatAstY),datLT,weig)
      cbr <- mvCumSum( as.vector( astY-hatAstY )*weig,datLT)/sqrt(n) 
      bootSamp <- rbind(bootSamp,cbind(Kn=max(abs(cbr)),Wn=mean(cbr^2)))
    }
  } else {
    objCall <- obj$call
    objmf <- getModelFrame(obj)
    for(i in 1:B)
    {
      objmf[,1] <- hatY + hatEps * rWildBoot(n)
      objCall$data <- objmf
      ##TODO: compIntRegProc(residuals(eval(objCall)),datLT,weig)
      cbr <- mvCumSum(as.vector(getResiduals(eval(objCall)))*weig,datLT)/sqrt(n)
      bootSamp <- rbind(bootSamp,c(Kn=max(abs(cbr)),Wn=mean(cbr^2)))
    }
  }
  return( bootSamp )
}


plotIntRegProc <- function(y,x,weig=rep(1,length(y)),ADD=FALSE,...)
##
#
{
  print(match.call())
  n <- length(y)
  stopifnot(n==length(x),n==length(weig))
  xord <- order(x)
  cumsum(y[xord]*weig[xord])/sqrt(n)
  if( ADD ) {
    points(x[xord],y[xord],...)
  } else {
    plot(c(x[1],x[xord]),c(0,y[xord]),type="n",main="",xlab="",ylab="")
    lines(range(x),c(0,0),col="black",lty=4) 
    points(x[xord],y[xord],...)
  }
}







#   getLessThan.R
#     Comparison needed to compute multivariate Integrated Residuals Tests.
#   implementaciones
#   NOTAS:
#     Mejorar con algoritomo de ordenación !!
#   ERRORES:


mvPartOrd <- function(x1,x2)
##OK
{
  stopifnot(is.vector(x1),is.vector(x2),length(x1)==length(x2))
  return( all(x1<=x2) )
}


getLessThan <- function(x,d)
##OK
#  x = a matrix
#  d = a matrix
#  RETURN
#    The list of elements in d that are less that x[i,]
#  for every i
{
  x <- as.matrix(x)
  d <- as.matrix(d)
  res <- rep(list(NULL),nrow(x))
  for(i in 1:nrow(d))
    for(j in 1:nrow(x))
      if( mvPartOrd(d[i,],x[j,]) )
        res[[j]] <- c(res[[j]],i)
  return(res)
}


mvCumSum <- function(x,ord)
##OK
#  x = data to sum up.
#  ord = order list, the result of 'getLessThan()'
{
  stopifnot(is.vector(x))
  res <- NULL
  for(iv in ord)
    if(is.null(iv))
      res <- c(res,0)
    else
      res <- c(res,sum(x[iv]))
  return( res )
}



# intRegGOF.R
#   Interface R to computed integrated regression.
# NOTAS:
#   - Does the order of variables play some role on the result ??
#   - TODO: check all covars are in data !!
# ERRORES:




.intRegGOFOptions <- list(B=999)


intRegGOF <- function(obj,covars=NULL,B=499,LINMOD=FALSE)
##
# obj = a lm, glm or nls fit object.
# covars = continuous variates uased to compute Int. Reg. 
# B = bootstrap sample size
# LINMOD = when true Linear Model matrix fitting equation is used.
##  OBS: nlm is not checked
{
  ##  
  irmc <- match.call()
  objCall <- obj$call
  ##  TODO: stopifnot if data argument is not set or check code !!
  objmf <- getModelFrame(obj)
  hatEps <- getResiduals(obj)
  hatY <- fitted(obj)
  weig <- getModelWeights(obj)
  if( is.null(covars) )
    covars <- getModelCovars(obj)
  else if( !is.vector(covars) )## assume it is a formula ~...
    covars <- all.vars(delete.response(covars))
## TODO: sólo las variables continuas !!!
  datXX <- eval( obj$call$data,envir=parent.frame() )
  covars <- getContVar(datXX,covars)
  if( length(covars)>0 ) { 
    datXX <- datXX[,covars]
    datLT <- getLessThan(as.matrix(datXX),as.matrix(datXX))
    bootSamp <- compBootSamp(obj,datLT,B,LINMOD)
    datIntErr <- mvCumSum( as.vector(hatEps)*weig,datLT)/sqrt(nrow(objmf)) 
    datKn <- max(abs(datIntErr))
    datWn <- mean(datIntErr^2)
    KnpVal <- mean(bootSamp[,"Kn"] > datKn)
    WnpVal <- mean(bootSamp[,"Wn"] > datWn)
  } else {
    # when y~1 we cannot do anything
    bootSamp <- datLT <- datIntErr <- datKn <- datWn <- KnpVal <- WnpVal <- NA
  }
  res <- list(  call=irmc,regObj=deparse(substitute(obj)),regModel=objCall,
                p.value=c(Kn=KnpVal,Wn=WnpVal), datStat=c(Kn=datKn,Wn=datWn),
                covars=covars,intErr=datIntErr,xLT=datLT,bootSamp=bootSamp)
  class(res) <- "intRegGOF"
  return( res )
}


print.intRegGOF <- function(x,...)
##
# x = an IntRegGOF object
{
  #cat("Integrated Regression Output\n")obj$irmc
  print(x$call)
  cat("Model Fit Call:\n  ")
  print( x$regModel )
  cat("Covariates: ",paste(x$covars,collapse=", "),".\n")
  rr <- cbind(value=x$datStat,p.value=x$p.value)
  rownames(rr) <- c("K","W^2")##c("max(abs(...))","mean(...^2)")
  print(rr)
}


plotAsIntRegGOF <- function(obj,covar=1,ADD=FALSE,...)
##
# obj = an lm, glm or nls object.
# covar = variable name or number whose Int. Reg. is computed. Number reference
#       vars in the model frame, while names refer to data
# add = if TRUE the plot is added to existing plot.
# ... = further parameters to plot
#   RETURN: Marginal Unidimensional plot for the Int. Residual function
{
  mc <- as.list(match.call())[-(1:4)] 
  if ( is.character(covar) ) {    ##TODO: Not well implemented
    stopifnot( length(covar)==1 )
    data <- eval( obj$call$data,envir=parent.frame()  )
    stopifnot( is.numeric(x <- data[,covar]) )
  } else if ( is.numeric(covar) && length(covar)==1 ) {
    stopifnot( covar>=0 )
    datXX <- getModelFrame(obj)
    covar <- names(datXX)[covar+1]
    x <- datXX[,covar]
  } else if ( is.numeric(covar) && length(covar)==length(getResiduals(obj)) ) {
    x <- covar
    covar <- deparse(substitute(covar))
  } else if ( class(covar)=="formula" ) {
    covar <- all.vars(covar)
    stopifnot( length(covar)==1 )
    data <- eval( obj$call$data,envir=parent.frame()  )
    stopifnot( is.numeric(x <- data[,covar]) )
  } else {
    ## PCA, bivariate stuff here ???
    stop(paste("'covar' has an improper value."))
  }
  if(ADD) {
    plotIntRegProc(getResiduals(obj),x,getModelWeights(obj),ADD=TRUE,...)
  } else {
    ## check for xlab, ylab, main
    print(mc)
    xlab <- if ( is.null(mc$xlab) ) covar else mc$xlab
    ylab <- if ( is.null(mc$ylab) ) "Int. Reg." else mc$ylab
    main <- if ( is.null(mc$main) ) obj$call else mc$main
    ## plotting instruction
    plotIntRegProc(getResiduals(obj),x,getModelWeights(obj),
                 ADD=FALSE,main=main,xlab=xlab,ylab=ylab,...)
  }
  invisible()
}


pointsAsIntRegGOF <- function(obj,covar=1,...)
{
  plotAsIntRegGOF(obj,covar=covar,ADD=TRUE,...)
  invisible()
}


linesAsIntRegGOF <- function(obj,covar=1,...)
{ 
  plotAsIntRegGOF(obj,covar=covar,ADD=TRUE,type="l",...)
  invisible()
}


anovarIntReg <- function(objH0,...,covars=NULL,B=499,
    LINMOD=FALSE,INCREMENTAL=FALSE)
##
# objH0 = a lm, glm or nls fit object. 
# ... = a number of lm, glm or nls fit objects.
# covars = continuous variates uased to compute Int. Reg. When NULL it is
#   obtained as the max. number of different covariates in all tested models 
# B = 
# LINMOD = when TRUE Linear Model matrix fitting equation is used.
##  choose between these two models from null distribution, the point 
##  of view of objH0 
##  make an interface that allows for multiple model comparison
##  TODO: INCREMENTAL should be coded
{
  nm <- 1
  if (is.null(covars))
  {
    covars <- getModelCovars(objH0)
    for(m in list(...)) {
      nm <- nm + 1
      cv <- getModelCovars(m)
      cvincov <- cv %in% covars 
      if ( !all(cvincov) )
        covars <- c(covars,cv[!cvincov])
    }
  } else if( !is.vector(covars) )## assume it is a formula ~...
    covars <- all.vars(delete.response(covars))
## TODO: sólo las variables continuas !!!
  # get weigths
  weig <- getModelWeights(objH0)
  # get data for cont. covars
  datXX <- eval(objH0$call$data,envir=parent.frame())
  datXX <- datXX[,getContVar(datXX,covars)]  
  # get LT struct in covars
  datLT <- getLessThan(as.matrix(datXX),as.matrix(datXX))
  # compute boot samp for objH0
  bs <- compBootSamp(objH0,datLT,B,LINMOD)
  # compute emp. proc for the rest of models
  res <- NULL
  for(m in list(objH0,...)) {
    ep <- compIntRegProc(resid(m),datLT,weig)
    epK <- max(abs(ep))
    epW <- mean(ep^2)
    KnpVal <- 2*min( mean( bs[,"Kn"] < epK ), mean( bs[,"Kn"] > epK ) )
    WnpVal <- 2*min( mean( bs[,"Wn"] < epW ), mean( bs[,"Wn"] > epW ) )
    if( INCREMENTAL && is.null(res) )
       bs <- compBootSamp(m,datLT,B,LINMOD)
    res <- rbind(res, c(epK,KnpVal,epW,WnpVal))
  }
  # compute table of stats and p.values
  mf <- lapply( list(objH0,...), function(x) x$call )
  tt <- "Integrated Regression Analysis of Variability Table:\n"
  if(INCREMENTAL) {
    tt <- paste(tt," Incremental Test Mode\n")
  } else {
    tt <- paste(tt," Reference Test Mode\n")
  }
  tn <- paste("  Model ", format((1:nm)-1), ": ", mf, 
              sep = "", collapse = "\n")
  tn <- paste(tn,"\nCovariables:",paste(covars,collapse=", "),"\n")
  colnames(res) <- c("K","P(>K)","W","P(>W)")
  rownames(res) <- paste("Model ", format((1:nm)-1), ": ", sep = "")
  return( structure(res, heading = c(tt, tn), class = c("anovarIntReg")) )
}



print.anovarIntReg <- function(x,...)
{
  cat( attr(x,"heading"), "\n", sep="" )
  printCoefmat(x,...)
}


# plot.anovarIntReg <- function(...)
# {
  # dibujar todos los procesos, con leyenda y tal 
# }


# wildBoot.R
#   utilidades para Wild Bootstrap.
# NOTAS:
# ERRORES:




## hMean <- function(x) return( 1/mean(1/x) )
## relVar <- function(x) return( mean( ((x-hMean(x))/x)^2 ) )


rWildBoot <- function(n)
##OK
#  Generate WILD WILD bootstrap distribution. 
#  X t.q. E[X]=0, E[X^2]=1, E[X^3]=1.
{
  B <- (1+sqrt(5))/2
  A <- (1-sqrt(5))/2
  p <- (5+sqrt(5))/10
  V <- ( A - B ) * rbinom(n,1,p) + B
  return( V )
}













