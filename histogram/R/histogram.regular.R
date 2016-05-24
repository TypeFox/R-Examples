`histogram.regular` <- function( y, penalty="br", breaks=NULL, control=list(), right=TRUE, verbose=TRUE, plot=TRUE, yvarname="y" ) {

  # check penalty-parameter 
  penalty = tolower( penalty )
  if ( penalty=="default" )
    penalty = "br"
  if ( ! any( penalty==c("br","aic","bic","cv","nml","sc","mdl")  ) ) {
    warning( "Penalty '", penalty, "' not supported for regular histograms - using 'br'", call.=FALSE )
    penalty = "br"
  }

  # check data vector
  if ( is.null(y) )
    stop( "no data in vector y" )

  lhvalue<-(-Inf)  
  y <- sort(y)
  n <- length(y)

  # control defaults 
  cont2 <- list( cvformula=1, p=1, g1=1, g2=1, g3=-1, mincount=0, maxbin=1000)
  if (penalty=="aic")
    cont2 <- list( cvformula=1, p=1, alpha=1, g1=1, g2=1, g3=-1, mincount=0, maxbin=1000)
  if (penalty=="bic")
    cont2 <- list( cvformula=1, p=1, alpha=0.5, g1=1, g2=1, g3=-1, mincount=0, maxbin=1000)

  # replace defaults with values set by the user
  cont2[names(control)]<-control
  
  # check g1/g2/g3
  if ( cont2$g1 <= 0 )
    stop( "g1 must be greater than 0" )
  if ( cont2$g2 < 0 )
    stop( "g2 must be greater than or equal to 0" )
  if ( cont2$g3 == -Inf )
    stop( "g3 must be a real number or Inf, but not -Inf" )


#check whether user specified breaks parameter
  if ( ! is.null(breaks) )
    Dmax<-floor(breaks)
  else  
    Dmax <- floor(cont2$g1*n^cont2$g2*log(n)^cont2$g3)

  # limit number of bins to maxbin
  Dmax = min( Dmax, cont2$maxbin ) 

  likelihood <- rep(0,Dmax)
  if ( verbose ) 
    message(paste("Building regular histogram with maximum number of bins ",Dmax,".",sep=""))

  # set mincount defaults
  if ( is.null(cont2$mincount) )
    cont2$mincount <- 0
  if ( ! any( names(control) == "mincount" ) ) {
    if ( penalty=="mdl" )
      cont2$mincount <- 1
    if ( (penalty=="cv")&&(cont2$cvformula==3) )
      cont2$mincount <- 2
  }
  # check mincount range
  cont2$mincount <- floor( cont2$mincount )
  if ( (penalty=="cv")&&(cont2$cvformula==3) ) 
    cont2$mincount <- max( cont2$mincount, 2 )
  if ( penalty=="mdl" ) 
    cont2$mincount <- max( cont2$mincount, 1 )  
  if ( cont2$mincount < 0 ) {
    warning( "mincount must not be negative - using 'mincount=0'", call.=FALSE )
    cont2$mincount <- 0
  }
  if ( cont2$mincount > floor(n/2) ) {
    warning( "mincount must not be greater than n/2 - using 'mincount=n/2'", call.=FALSE )
    cont2$mincount <- floor(n/2)
  }
    

  # check cvformula
  if ( !is.null(cont2$cvformula) )
    if ( cont2$cvformula!=1 && cont2$cvformula!=2 && cont2$cvformula!=3 )
      stop( "cvformula must be 1, 2 or 3" )
      
  # check p
  if ( !is.null(cont2$p) ) {
    if ( (cont2$p%%1!=0) || (cont2$p<1) || (cont2$p>n-1) )
      stop( "p must be an integer between 1 and (n-1)" )
    if ( cont2$p!=1 && cont2$cvformula!=2 ) {
      # if p!=1: always use cvformula 2 
      warning( "For p!=1 cvformula must be set to 2 - using 'cvformula=2'", call.=FALSE )
      cont2$cvformula = 2
    }
  }

  
  ## maximum Likelihood options
  is.ml<-FALSE
  if ( penalty=="br" ) {
    pen <- (1:Dmax)+(log(1:Dmax))^2.5
    is.ml <- TRUE
  }
  if ( penalty=="aic" )  {
    pen<-cont2$alpha*(1:Dmax)
    is.ml<-TRUE
  }
  if ( penalty=="bic" ) {
    pen<-cont2$alpha*log(n)*(1:Dmax)
    is.ml<-TRUE
  }
  if ( penalty=="nml" )  {
    m <- 1:Dmax
    bt <- gamma(1/2)/beta( (m-1)/2, 1/2 )
    pen <- ( (m-1)/2 *log(n/2) 
            + log( sqrt(pi) ) - lgamma(m/2)
            + (sqrt(2)*m)/(3*sqrt(n)) * bt
            + 1/n * ((3+m*(m-2)*(2*m+1))/36)
            - m^2/(9*n)* bt^2
           )
    is.ml<-TRUE
  }

  if ( is.ml==TRUE ) {
	if (verbose) message("- Choosing number of bins via maximum likelihood with ",toupper(penalty)," penalty.",sep="") 
   for ( D in 1:Dmax ) {
      reghist <- hist( y, breaks=(y[1] + ((0:(D))/(D))*(y[n]-y[1])), right=right, plot=FALSE )
#      reghist <- hist( y, breaks=seq(y[1],y[n], length.out=(D+1)), right=right, plot=FALSE )
      like <- rep( 0,D )
      like2 <- log( reghist$density )
      like[is.finite(like2)] <- like2[is.finite(like2)]
      likelihood[D] <- sum(reghist$counts*like)
	if (min(reghist$counts)<cont2$mincount) likelihood[D]<--Inf
   }

    penlike <- likelihood-pen
    Dopt <- which.max(penlike)
    lhvalue<-max(penlike)+1
  }
  else {
	if (penalty=="cv") {
    for ( D in 1:Dmax ) {
      reghist <- hist( y, breaks=(y[1] + ((0:(D))/(D))*(y[n]-y[1])), right=right, plot=FALSE )
#      reghist <- hist( y, breaks=seq(y[1],y[n], length.out=(D+1)), right=right, plot=FALSE )
      if ( cont2$cvformula==2 )
        likelihood[D] <- D*(n-cont2$p+1)/n*sum(reghist$counts^2)-(2*n-cont2$p)*D
      else
        if ( cont2$cvformula==3 ) {
	if (min(reghist$counts)<cont2$mincount) likelihood[D]<--Inf
          else likelihood[D] <- sum(reghist$counts*log(reghist$counts-1))+n*log(D)
        }
        else
          likelihood[D]<-D*(n+1)/n^2*sum(reghist$counts^2)-2*D
	if (min(reghist$counts)<cont2$mincount) likelihood[D]<--Inf
    }
	if (verbose) message(paste("- Choosing number of bins via leave-",cont2$p,"-out cross validation. Using formula ",cont2$cvformula,".",sep=""))
    Dopt<-which.max(likelihood)
	}
	if (penalty=="sc") {
		if (verbose) message("- Choosing number of bins via stochastic complexity (SC) criterion.",sep="") 
	   for ( D in 1:Dmax ) {
	      reghist <- hist( y, breaks=(y[1] + ((0:(D))/(D))*(y[n]-y[1])), right=right, plot=FALSE )
#	      reghist <- hist( y, breaks=seq(y[1],y[n], length.out=(D+1)), right=right, plot=FALSE )
     	      likelihood[D] <- sum(lfactorial(reghist$counts))-lchoose(D+n-1,D-1)+n*log(D)
		if (min(reghist$counts)<cont2$mincount) likelihood[D]<--Inf
	    }
	    Dopt <- which.max(likelihood)
	}
	if (penalty=="mdl") {
		if (verbose) message("- Choosing number of bins via minimum description length (MDL) criterion.") 
	   for ( D in 1:Dmax ) {
	      reghist <- hist( y, breaks=(y[1] + ((0:(D))/(D))*(y[n]-y[1])), right=right, plot=FALSE )	   
#	      reghist <- hist( y, breaks=seq(y[1],y[n], length.out=(D+1)), right=right, plot=FALSE )

		if (min(reghist$counts)<cont2$mincount) likelihood[D]<--Inf
		else likelihood[D]<-sum((reghist$counts-0.5)*log(reghist$counts-0.5))-(n-D/2)*log(n-D/2)+n*log(D)-0.5*D*log(n)
	    }
	    Dopt <- which.max(likelihood)
	}


  }
  

  # create histogram

	if (verbose) message(paste("- Number of bins chosen: ",Dopt,".\n\n",sep=""))
  H <- hist( y, breaks=(y[1] + ((0:(Dopt))/(Dopt))*(y[n]-y[1])), right=right, plot=FALSE )
#  H <- hist( y, breaks=seq(y[1], y[n], length.out=(Dopt+1)), right=right, plot=FALSE )

  # Bugfix: Name of y-var gets lost above - reset it.
  H$xname = yvarname


  if ( plot )
    plot( H, freq=FALSE ) 

  return(list(H=H,lhvalue=lhvalue) )
}

