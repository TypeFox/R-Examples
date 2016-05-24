plotPostPred <-
function(BESTobj, nCurvesToPlot = 30) {
  # This function plots the posterior predictive distribution and the data. 
  # Description of arguments:
  # BESTobj is mcmc.list object of the type returned by function BESTmcmc.

  # Sanity checks:
  if(!inherits(BESTobj, "data.frame"))
    stop("BESTobj is not a valid BEST object")
  if(ncol(BESTobj) == 3 && all(colnames(BESTobj) == c("mu","nu","sigma"))) {
    oneGrp <- TRUE
  } else if (ncol(BESTobj) == 5 && all(colnames(BESTobj) == c("mu1", "mu2","nu","sigma1","sigma2"))) {
    oneGrp <- FALSE
  } else {
    stop("BESTobj is not a valid BEST object")
  }


  # mcmcChain <- as.matrix(BESTobj)
  y1 <- attr(BESTobj, "data")$y1
  y2 <- attr(BESTobj, "data")$y2

  # Set up window and layout:
  oldpar <- par(mar=c(3.5,3.5,2.5,0.5), mgp=c(2.25,0.7,0), "mfrow") 
    on.exit(par(oldpar))
  if(!oneGrp)
    par(mfrow=2:1)
 
  # Select thinned steps in chain for plotting of posterior predictive curves:
  stepIdxVec <- seq(1, NROW( BESTobj ), length= nCurvesToPlot)
  toPlot <- BESTobj[stepIdxVec, ]

  if(oneGrp)  {
    mu1 <- toPlot$mu
    sigma1 <- toPlot$sigma
    y2 <- mu2 <- sigma2 <- NULL
  } else {
    mu1 <- toPlot$mu1
    mu2 <- toPlot$mu2
    sigma1 <- toPlot$sigma1
    sigma2 <- toPlot$sigma2
  }
  nu <- toPlot$nu

  if(is.null(y1) && is.null(y2)) {
    xLim <- range(mu1, mu2)
  } else {
    xRange <- range(y1, y2)
    xLim <- c( xRange[1]-0.1*diff(xRange) , 
              xRange[2]+0.1*diff(xRange) )
  }
  xVec <- seq(xLim[1], xLim[2], length=200)
  maxY <- max( dt(0, df=max(nu)) / min(sigma1, sigma2) )

  # Plot data and smattering of posterior predictive curves:
  isErr <- try(plotDataPPC(y=y1, mu=mu1, sigma=sigma1, nu=nu, xVec=xVec, maxY=maxY), silent=TRUE)
  if(inherits(isErr, "try-error"))
    text(min(xVec), maxY, "Cannot plot data here.", pos=4, col='red')
  if(oneGrp) {
    title(main="Data w. Post. Pred.")
    if(!is.null(y1))
      text( max(xVec) , maxY , bquote(N ==.(length(y1))) , adj=c(1.1,1.1) )
  } else {
    title(main="Data Group 1 w. Post. Pred.")
    if(!is.null(y1))
      text( max(xVec) , maxY , bquote(N[1]==.(length(y1))) , adj=c(1.1,1.1) )
  }
  if(!oneGrp) {
    isErr <- try(plotDataPPC(y=y2, mu=mu2, sigma=sigma2, nu=nu, xVec=xVec, maxY=maxY), silent=TRUE)
    if(inherits(isErr, "try-error"))
      text(min(xVec), maxY, "Cannot plot data here.", pos=4, col='red')
   title(main="Data Group 2 w. Post. Pred.")
    if(!is.null(y2))
      text( max(xVec) , maxY , bquote(N[2]==.(length(y2))) , adj=c(1.1,1.1) )
  }

  return(invisible(NULL))
}
