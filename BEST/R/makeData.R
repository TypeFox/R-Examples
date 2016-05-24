makeData <-
function( mu1 , sd1 , mu2=NULL , sd2=NULL , nPerGrp ,
                     pcntOut=0 , sdOutMult=2.0 ,
                     rnd.seed=NULL, showPlot=TRUE ) {
  # Auxilliary function for generating random values from a
  # mixture of normal distibutions.

  oneGrp <- is.null(mu2) || is.null(sd2)
  if(!is.null(rnd.seed)){set.seed(rnd.seed)} # Set seed for random values.
  nOut = ceiling(nPerGrp*pcntOut/100)        # Number of outliers.
  nIn = nPerGrp - nOut                       # Number from main distribution.

  # Sanity checks
  if ( pcntOut > 100 || pcntOut < 0 )
    stop("pcntOut must be between 0 and 100.")
  if ( pcntOut > 0 && pcntOut < 1 )
    warning("pcntOut is specified as percentage 0-100, not proportion 0-1.")
  if ( pcntOut > 50 )
    warning("pcntOut indicates more than 50% outliers; did you intend this?")
  if ( nOut < 2 && pcntOut > 0 )
    stop("Combination of nPerGrp and pcntOut yields too few outliers.")
  if ( nIn < 2  )
    stop("Combination of nPerGrp and pcntOut yields too few non-outliers.")

  sdN = function( x )
    sqrt( mean( (x-mean(x))^2 ) )
  exactify <- function(x, mean, sd) # Scale to exact mean and sdN.
    (x - mean(x))/sdN(x) * sd + mean

  y1 <- exactify(rnorm(n=nIn), mu1, sd1)  # main distribution
  if(nOut > 0)
    y1 <- c(y1, exactify(rnorm(n=nOut), mu1, sd1*sdOutMult)) # outliers
  if(oneGrp) {
    y2 <- NULL
  } else {
    y2 <- exactify(rnorm(n=nIn), mu2, sd2)
    if(nOut > 0)
      y2 <- c(y2, exactify(rnorm(n=nOut), mu2, sd2*sdOutMult))
  }
  #
  if(showPlot) {
    # Set up window and layout:
    opar <- par(mfrow=c(1,1)) ; on.exit(par(opar))
    if(!oneGrp)
      par(mfrow=2:1)
      # layout(matrix(1:2,nrow=2))
    histInfo = hist( y1 , main="Simulated Data" , col="pink2" , border="white" ,
                     xlim=range(c(y1,y2)) , breaks=30 , prob=TRUE )
    text( max(c(y1,y2)) , max(histInfo$density) ,
          bquote(N==.(nPerGrp)) , adj=c(1,1) )
    xlow=min(histInfo$breaks)
    xhi=max(histInfo$breaks)
    xcomb=seq(xlow,xhi,length=1001)
    lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1)*nIn/(nIn+nOut) +
      dnorm(xcomb,mean=mu1,sd=sd1*sdOutMult)*nOut/(nIn+nOut) , lwd=3 )
    lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1)*nIn/(nIn+nOut) ,
           lty="dashed" , col="blue", lwd=3)
    lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1*sdOutMult)*nOut/(nIn+nOut) ,
           lty="dashed" , col="red", lwd=3)
    if(!oneGrp)  {
      histInfo = hist( y2 , main="" , col="pink2" , border="white" ,
                       xlim=range(c(y1,y2)) , breaks=30 , prob=TRUE )
      text( max(c(y1,y2)) , max(histInfo$density) ,
            bquote(N==.(nPerGrp)) , adj=c(1,1) )
      xlow=min(histInfo$breaks)
      xhi=max(histInfo$breaks)
      xcomb=seq(xlow,xhi,length=1001)
      lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2)*nIn/(nIn+nOut) +
        dnorm(xcomb,mean=mu2,sd=sd2*sdOutMult)*nOut/(nIn+nOut) , lwd=3)
      lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2)*nIn/(nIn+nOut) ,
             lty="dashed" , col="blue", lwd=3)
      lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2*sdOutMult)*nOut/(nIn+nOut) ,
             lty="dashed" , col="red", lwd=3)
    }
  }
  #
  return( list( y1=y1 , y2=y2 ) )
}
