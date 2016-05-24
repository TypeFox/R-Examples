plot.lmWinsor <- function(x, n=101, lty=1:9, col=1:9,
         lwd=c(2:4, rep(3, 6)), lty.y=c('dotted', 'dashed'),
         lty.x = lty.y, col.y=1:9, col.x= col.y, lwd.y = c(1.2, 1),
         lwd.x=lwd.y, ...){
##
## 1.  One fit or several?
##
  clo <- class(x)
  {
    if(length(clo)<2){
      objNames <- names(x)
      objs <- x[objNames != 'call']
      obj1 <- objs[[1]]
    }
    else{
      objs <- list(x) 
      obj1 <- x     
    }
  }
  nfits <- length(objs) 
##
## 2.  Number of explanatory variables?
##
  mdly <- mdlx <- mdl <- formula(obj1)
  mdly[[3]] <- NULL
  mdlx[[2]] <- NULL
  xNames <- all.vars(mdlx)
  yName <- all.vars(mdly)
#
  kx <- length(xNames)
##
## 3.  Plot the data
##
  Xy <- model.frame(obj1)
  yi <- Xy[, yName] 
  {
    if(kx == 1){
      xx <- Xy[, xNames]
      xlab <- xNames
    }
    else{
      xx <- predict(obj1)
      xlab <- 'predicted' 
    }
  }
#
  dots <- list(...)
  if(!('xlab' %in% names(dots)))
     dots$xlab <- xlab
  if(!('ylab' %in% names(dots)))
    dots$ylab <- yName
  {
    if('xlim' %in% names(dots))
      xlim <- dots$xlim
    else {
      xlim <- range(xx)
      dots$xlim <- xlim
    }
  }  
#
  dots$x <- xx
  dots$y <- yi
  if(!('main' %in% names(dots))){
    cl <- x$call
    main <- paste(as.character(cl$formula)[c(2,1,3)], collapse=" ")
    if(class(cl$data)=='name')
      main <- paste(as.character(cl$data), main, sep=':  ') 
#
    if('trim' %in% names(cl)){
      trim <- paste(eval(cl$trim), collapse=", ")
      main <- paste(main, trim, sep=";  trim = ") 
    }
    dots$main <- main
  }
  do.call('plot', dots)
##
## 4.  Plot one line for each fit
##
#  4.1.  x for lines   
  if(kx==1){
    x. <- seq(xlim[1], xlim[2], length=n)
    dfx <- data.frame(x.)
    names(dfx) <- xNames 
  }
#  4.2.  Graphics parameters for the lines
  lty <- rep(lty, length=nfits) 
  col <- rep(col, length=nfits) 
  lwd <- rep(lwd, length=nfits)
#  
  lty.y <- rep(lty.y, length=nfits) 
  col.y <- rep(col.y, length=nfits) 
  lwd.y <- rep(lwd.y, length=nfits) 
#
  lty.x <- rep(lty.x, length=nfits) 
  col.x <- rep(col.x, length=nfits) 
  lwd.x <- rep(lwd.x, length=nfits)
#  4.3.  one line for each fit
    
  for(ifit in 1:nfits){
    Li <- objs[[ifit]]$lower
    Ui <- objs[[ifit]]$upper
    abline(h=c(Li[yName], Ui[yName]), lty=lty.y[ifit],
           col=col.y[ifit], lwd=lwd.y[ifit]) 
    if(kx==1){
      obji <- objs[[ifit]]
#      envi <- new.env()
#      assign('obj.', obji, envir=envi) 
#      assign('newdat.', dfx, envir=envi)
#      predList <- list(object=quote(obj.), newdata=quote(newdat.) )
      predList <- list(object=obji, newdata=dfx )
#
      y. <- do.call('predict', predList)
#      y. <- do.call('predict', predList, envir=envi)
#      
      lines(x., y., lty=lty[ifit], col=col[ifit],
            lwd=lwd[ifit])
      abline(v=c(Li[xNames], Ui[xNames]), lty=lty.x[ifit],
             col=col.x[ifit], lwd=lwd.x[ifit])
    }
    else{
      yp <- predict(objs[[ifit]])
      yo <- order(yi, yp)
      lines(yp[yo], yi[yo], lty=lty[ifit], col=col[ifit],
            lwd=lwd[ifit] )
      abline(v=c(Li[yName], Ui[yName]), lty=lty.x[ifit],
             col=col.x[ifit], lwd=lwd.x[ifit])
    }
  }
##
## 5.  Done
##
  invisible(NULL)
}
