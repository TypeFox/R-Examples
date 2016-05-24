plotAll <-
function(BESTobj, credMass=0.95,
                    ROPEm=NULL, ROPEsd=NULL, ROPEeff=NULL,
                    compValm=0, compValsd=NULL, compValeff=0,
                    showCurve=FALSE, ...) {
  # This function plots the posterior distribution (and data). It does not
  #   produce a scatterplot matrix; use pairs(...) for that.
  # Description of arguments:
  # BESTobj is BEST object of the type returned by function BESTmcmc.
  # ROPEm is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of means.
  # ROPEsd is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of standard deviations.
  # ROPEeff is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the effect size.
  # showCurve if TRUE the posterior should be displayed as a fitted density curve
  #   instead of a histogram (default).

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

  y1 <- attr(BESTobj, "data")$y1
  y2 <- attr(BESTobj, "data")$y2

  # Select thinned steps in chain for plotting of posterior predictive curves:
  chainLength <- NROW( BESTobj )
  nCurvesToPlot <- 30
  # stepIdxVec <- seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
  stepIdxVec <- seq(1, chainLength, length.out=nCurvesToPlot)

  if(oneGrp)  {
    mu1 <- BESTobj$mu
    sigma1 <- BESTobj$sigma
    y2 <- mu2 <- sigma2 <- NULL
  } else {
    mu1 <- BESTobj$mu1
    mu2 <- BESTobj$mu2
    sigma1 <- BESTobj$sigma1
    sigma2 <- BESTobj$sigma2
  }
  nu <- BESTobj$nu

  # Set up window and layout:
  # windows(width=6.0,height=8.0) # Better to use default plot window.
  oldpar <- par(mar=c(3.5,3.5,2.5,0.5), mgp=c(2.25,0.7,0), "mfrow") 
    on.exit(par(oldpar))
  if(oneGrp) {
    layout( matrix( c(3,3,4,4,5,5, 1,1,1,1,2,2) , nrow=6, ncol=2 , byrow=FALSE ) )
  } else {
    layout( matrix( c(4,5,7,8,3,1,2,6,9,10) , nrow=5, byrow=FALSE ) )
  }

  if(is.null(y1) && is.null(y2)) {
    xLim <- range(mu1, mu2)
  } else {
    xRange <- range(y1, y2)
    xLim <- c( xRange[1]-0.1*diff(xRange) , 
              xRange[2]+0.1*diff(xRange) )
  }
  xVec <- seq( xLim[1] , xLim[2] , length=200 )
  maxY <- max( dt( 0 , df=max(nu) ) /
    min(sigma1,sigma2) )


  # Plot data and smattering of posterior predictive curves:
  plotDataPPC(y=y1, mu=mu1[stepIdxVec], sigma=sigma1[stepIdxVec], nu=nu[stepIdxVec],
    xVec=xVec, maxY=maxY)
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
    plotDataPPC(y=y2, mu=mu2[stepIdxVec], sigma=sigma2[stepIdxVec], nu=nu[stepIdxVec],
      xVec=xVec, maxY=maxY)
    title(main="Data Group 2 w. Post. Pred.")
    if(!is.null(y2))
      text( max(xVec) , maxY , bquote(N[2]==.(length(y2))) , adj=c(1.1,1.1) )
  }

  # Plot posterior distribution of parameter nu:
  plotPost( log10(nu) , col="skyblue" ,  credMass=credMass,
                       showCurve=showCurve ,
                  xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , showMode=TRUE ,
                  main="Normality" ) #  (<0.7 suggests kurtosis)

  # Plot posterior distribution of parameters mu1, mu2, and their difference:
  xlim <- range( c( mu1 , mu2 ) )
  if(oneGrp) {
    plotPost( mu1 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                       showCurve=showCurve , ROPE=ROPEm, compVal=compValm,
                  xlab=bquote(mu) , main="Mean" , 
                  col="skyblue" )
  } else {
    plotPost( mu1 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                       showCurve=showCurve ,
                  xlab=bquote(mu[1]) , main=paste("Group",1,"Mean") , 
                  col="skyblue" )
    plotPost( mu2 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                         showCurve=showCurve ,
                    xlab=bquote(mu[2]) , main=paste("Group",2,"Mean") , 
                    col="skyblue" )
    plotPost( mu1-mu2 , compVal=0 ,  showCurve=showCurve , credMass=credMass,
                    xlab=bquote(mu[1] - mu[2]) , cex.lab = 1.75 , ROPE=ROPEm ,
                    main="Difference of Means" , col="skyblue" )
  }

  # Plot posterior distribution of param's sigma1, sigma2, and their difference:
  xlim <- range( c( sigma1 , sigma2 ) )
  if(oneGrp) {
    plotPost(sigma1, xlim=xlim, cex.lab = 1.75, credMass=credMass,
                       showCurve=showCurve, ROPE=ROPEsd, compVal=compValsd,
                  xlab=bquote(sigma) , main="Std. Dev." , 
                  col="skyblue" , showMode=TRUE )
  } else {
    plotPost( sigma1 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                       showCurve=showCurve ,
                  xlab=bquote(sigma[1]) , main=paste("Group",1,"Std. Dev.") , 
                  col="skyblue" , showMode=TRUE )
    plotPost( sigma2 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                         showCurve=showCurve ,
                    xlab=bquote(sigma[2]) , main=paste("Group",2,"Std. Dev.") , 
                    col="skyblue" , showMode=TRUE )
    plotPost( sigma1-sigma2 ,  credMass=credMass,
                         compVal=compValsd ,  showCurve=showCurve ,
                         xlab=bquote(sigma[1] - sigma[2]) , cex.lab = 1.75 , 
                         ROPE=ROPEsd ,
                 main="Difference of Std. Dev.s" , col="skyblue" , showMode=TRUE )
  }

  # Effect size for 1 group:
  if(oneGrp)  {
  effectSize = ( mu1 - compValm ) / sigma1
  plotPost( effectSize , compVal=compValeff ,  ROPE=ROPEeff , credMass=credMass,
                 showCurve=showCurve ,
                 xlab=bquote( (mu-.(compValm)) / sigma ),
                 showMode=TRUE , cex.lab=1.75 , main="Effect Size" , col="skyblue" )
  } else {
    # Plot of estimated effect size. Effect size is d-sub-a from 
    # Macmillan & Creelman, 1991; Simpson & Fitter, 1973; Swets, 1986a, 1986b.
    effectSize <- ( mu1 - mu2 ) / sqrt( ( sigma1^2 + sigma2^2 ) / 2 )
    plotPost( effectSize , compVal=compValeff ,  ROPE=ROPEeff , credMass=credMass,
                          showCurve=showCurve ,
                    xlab=bquote( (mu[1]-mu[2])
                      /sqrt((sigma[1]^2 +sigma[2]^2 )/2 ) ),
                showMode=TRUE , cex.lab=1.0 , main="Effect Size" , col="skyblue" )
  }
  # Or use sample-size weighted version:
  # Hedges 1981; Wetzels, Raaijmakers, Jakab & Wagenmakers 2009.
  # N1 = length(y1)
  # N2 = length(y2)
  # effectSize = ( mu1 - mu2 ) / sqrt( ( sigma1^2 *(N1-1) + sigma2^2 *(N2-1) )
  #                                    / (N1+N2-2) )
  # Be sure also to change BESTsummary function, above.
  # histInfo = plotPost( effectSize , compVal=0 ,  ROPE=ROPEeff ,
  #          showCurve=showCurve ,
  #          xlab=bquote( (mu[1]-mu[2])
  #          /sqrt((sigma[1]^2 *(N[1]-1)+sigma[2]^2 *(N[2]-1))/(N[1]+N[2]-2)) ),
  #          showMode=TRUE , cex.lab=1.0 , main="Effect Size" , col="skyblue" )
  return(invisible(NULL))
}
