`histogram.irregular` <- function( y, grid="data", breaks=NULL, penalty="penB", greedy=TRUE, right=TRUE, control=list(), verbose=TRUE, plot=TRUE, yvarname="y" ) {

  epsilon <- 1e-7

  penalty = tolower( penalty )
  grid = tolower(grid)

  # check penalty parameter
  if ( penalty=="default" )
    penalty = "penb"
  if ( ! any( penalty==c("penr","pena","penb","aic","bic","cv")  ) ) {
    warning( "Penalty '", penalty, "' not supported for irregular histograms - using 'penR'", call.=FALSE )
    penalty = "penr"
  }

  # check data vector
  if ( is.null(y) )
    stop( "no data in vector y" )
  lhvalue<-(-Inf)
  n <- length(y)
  a <- min(y)
  b <- max(y)
  y <- sort((y-a)/(b-a)) # normalize range to [0,1]

  if (verbose) 
    message("Building irregular histogram.")

  # control defaults 
  if ( penalty=="penr" )
    cont2 <- list( g1=1, g2=1, alpha=0.5, c=1 )
  if (penalty=="penb")
    cont2 <- list( g1=1, g2=1, alpha=1,   c=1 )
  if (penalty=="pena")
    cont2 <- list( g1=1, g2=1, alpha=0.5,   c=1,  k=2 )
  if (penalty=="aic")
    cont2 <- list( g1=1, g2=1, alpha=1 )
  if (penalty=="bic")
    cont2 <- list( g1=1, g2=1, alpha=0.5 )
  if (penalty=="cv")
    cont2 <- list( g1=1, g2=1, cvformula=1, p=1 )

  cont2$quanttype=7
  
  cont2$maxbin=1000
  if ( greedy == TRUE )
    cont2$maxbin=Inf

	if (grid=="data") 
	  cont2$g3=Inf
	if (grid=="regular") 
	  cont2$g3=-1
	if (grid=="quantiles") 
	  cont2$g3=-1

  cont2$between<-FALSE

  # replace control defaults with values set by the user
  cont2[names(control)]<-control
  
  # check g1/g2/g3
  if ( cont2$g1 <= 0 )
    stop( "g1 must be greater than 0" )
  if ( cont2$g2 < 0 )
    stop( "g2 must be greater than or equal to 0" )
  if ( cont2$g3 == -Inf )
    stop( "g3 must be a real number or Inf, but not -Inf" )  

### if user supplies breaks argument, use this for G
### breakes argument now checked in the wrapper function

  if ( !is.null(breaks) )
    G <- floor( breaks )
  else
    G <- cont2$g1*n^cont2$g2*log(n)^cont2$g3


  # check controls
  # alpha
  if ( !is.null(cont2$alpha) )
    if ( cont2$alpha < 0 )
      stop( "alpha must be greater than or equal to 0" )
  # c
  if ( !is.null(cont2$c) )
    if ( cont2$c < 0 )
      stop( "c must be greater than or equal to 0" )
  # k
  if ( !is.null(cont2$k) )
    if ( cont2$k < 1 )
      stop( "k must be greater than or equal to 1" )
  # cvformula
  if ( !is.null(cont2$cvformula) )
    if ( cont2$cvformula!=1 && cont2$cvformula!=2)
      stop( "cvformula must be 1 or 2" )
  # p
  if ( !is.null(cont2$p) ) {
    if ( (cont2$p%%1!=0) || (cont2$p<1) || (cont2$p>n-1) )
      stop( "p must be an integer between 1 and (n-1)" )
    if ( cont2$p!=1 && cont2$cvformula!=2 ) {
      # falls p!=1: cvformula 2 verwenden
      warning( "For p!=1 'cvformula=2' is used", call.=FALSE )
      cont2$cvformula = 2
    }
  }
  # mincount                  
  mincount=0
  if ( (penalty=="cv")&&(cont2$cvformula==3) )
    mincount <- 2
  
  # check grid parameter
  if ( ! any( grid==c("data","regular","quantiles")  ) ) {
    warning( "Grid '", grid, "' not supported - using 'data'", call.=FALSE )
    grid = "data"
  }
  
  if ( grid=="regular" ) {
    if ( length(breaks)==1 )
       BL <- (0:((breaks+1)-1)) / (((breaks+1)-1))
#      BL<-seq(0,1,length.out=(breaks+1))
    if ( is.null(breaks) )
       BL <- (0:(max(2,floor(G)+1)-1)) / (max(2,floor(G)+1)-1)
#      BL<-seq(0,1,length.out=max(2,floor(G)+1))
    if (verbose) 
      message(paste("- Using a regular grid with ",length(BL)-1," bins as finest grid.",sep=""))
  }

  if (grid=="quantiles") {
    if (length(breaks)==1)
       bk <- (0:((breaks+1)-1)) / (((breaks+1)-1))
#      bk<-seq(0,1,length.out=(breaks+1))
    if (is.null(breaks))
       bk <- (0:(max(2,floor(G)+1)-1)) / (max(2,floor(G)+1)-1)
#      bk<-seq(0,1,length.out=max(2,floor(G)+1))
    BL<-c(0,quantile(y,probs=bk[-c(1,length(bk))],type=cont2$quanttype),1)
	  BL<-unique(sort(BL))
    if (verbose) 
      message("- Using a regular quantile grid with ",length(BL)-1," bins as finest grid.",sep="")
  }

  if (grid=="data") {
    if (cont2$between) BL<-c(y[1]-epsilon,y[-n]+diff(y)/2,y[n]+epsilon)
    else {
      if (!right)
        BL <- c(y[1:n-1]-epsilon,(y[n-1]+y[n])/2,y[n]+epsilon)
      else
        BL <- c(y[1]-epsilon,(y[1]+y[2])/2,y[2:n]+epsilon)
  	}
  	BL<-unique(sort(BL))
    if (verbose) 
      message("- Using finest grid based on observations.")
  }

  # specify weightfunctions for dynamical programming part...
  
  if (penalty=="penr") {
	if (verbose) message("- Choosing number of bins via maximum likelihood with ",toupper(penalty)," penalty.",sep="") 
    mini<-FALSE
    weightfunction <- function(i,j) {
      dN = NL[j]-NL[i]
      dBL = BL[j]-BL[i]
      # compute the loglikelihood part of interval [B(i),B(j)[
      if (dN==0)
        W<-0
      else
        if ((dBL<1/G)&(grid=="data"))
          return(-Inf) # remove small bins
        else
          W <- dN*log(dN/dBL)
       W <- W-cont2$alpha*dN/dBL/n  ## use of random penalty
	 if (dN<mincount) W<- -Inf
       return(W)
    }
  }

  if ( penalty=="pena" || penalty=="penb" || penalty=="aic" || penalty=="bic" ) {
	if (verbose) message("- Choosing number of bins via maximum likelihood with ",toupper(penalty)," penalty.",sep="") 
    mini<-FALSE
    weightfunction<-function(i,j) {
      dN = NL[j]-NL[i]
      dBL = BL[j]-BL[i]
      # compute the loglikelihood part of interval [B(i),B(j)[
      if (dN==0)
        W<-0
      else
        if ((dBL<1/G)&(grid=="data"))
          return(-Inf) # remove small bins
        else
          W <- dN*log(dN/dBL)
	 if (dN<mincount) W<- -Inf
      return(W)
    }
  }

  if (penalty=="cv") {
	if (verbose) message(paste("- Choosing number of bins via leave-",cont2$p,"-out cross validation. Using formula ",cont2$cvformula,".",sep=""))
    mini<-TRUE
    if ( cont2$cvformula==2 )
      weightfunction<-function(i,j) {
        dN = NL[j]-NL[i]
        dBL = BL[j]-BL[i]
        # compute the cv part for interval
        if (dN==0)
          W<-0
        else
          if ((dBL<1/G)&(grid=="data"))
            return(Inf) # remove small bins
          else
            W <- ((2*n-cont2$p)*dN/(n*dBL)-(n-cont2$p+1)*dN^2/(dBL*n) )/((n-1)*(n-cont2$p))
	 if (dN<mincount) W<- Inf
        return(W)
      }
    else
      weightfunction<-function(i,j) {
        dN = NL[j]-NL[i]
        dBL = BL[j]-BL[i]
        # compute the cv part for interval
        if (dN==0)
          W<-0
        else
          if ((dBL<1/G)&(grid=="data"))
            return(Inf) # remove small bins
          else
            W <- (2*dN/(n*dBL)-(n+1)*dN^2/(dBL*n^2) )/(n-1)
	 if (dN<mincount) W<- Inf
       return(W)
      }
  }


  NL<-calcNL(y,BL,right)


  # Greedy
  if ( greedy==TRUE ) {
	D <- length(BL)-1 
      binmax<-ceiling(min(max(100,D^(1/3)),D))
	if ((D>=2)&(D>binmax)){
    if (verbose)
      message(paste("- Using greedy procedure to recursively build a finest partition with at most ",binmax," bins.",sep=""))
    BL<-histgreedy(BL,NL,n,binmax,verbose=verbose)
    NL<-calcNL(y,BL,right)
	}
  }

  D <- length(BL)-1
  #D <- min(length(NL)-1,D)
  #bounds<-BL


  # optimize & penalize
  if (penalty=="penr") {
    pen <- cont2$c*logchoose(n-1,D-1) + log(1:D)^2.5
    # compute dynamically the max likelihood for irregular histograms with
    # d bins, d varying between 1 and D
    tmp <- DynamicExtreme(weightfunction,n=length(BL)-1,D=D,mini=FALSE,msg=verbose)
    D0 <- which.max(tmp$extreme-pen)
    lhvalue<-max(tmp$extreme-pen)-n*log(b-a)-n*log(n)+cont2$alpha
  }
  if (penalty=="penb") {
    pen <- cont2$c*logchoose(n-1,D-1) + cont2$alpha*(0:(D-1)) + log(1:D)^2.5
    # compute dynamically the max likelihood for irregular histograms with
    # d bins, d varying between 1 and D
    tmp <- DynamicExtreme(weightfunction,n=length(BL)-1,D=D,mini=FALSE,msg=verbose)
    D0 <- which.max(tmp$extreme-pen)
    lhvalue<-max(tmp$extreme-pen)-n*log(b-a)-n*log(n)
  }
  if (penalty=="pena") {
    pen <- cont2$c*logchoose(n-1,D-1) + cont2$alpha*(0:(D-1)) + cont2$c*cont2$k*log(1:D)+
    2*sqrt(cont2$c*cont2$alpha*(0:(D-1))*(logchoose(n-1,D-1)+cont2$k*log(1:D)) )
    # compute dynamically the max likelihood for irregular histograms with
    # d bins, d varying between 1 and D
    tmp <- DynamicExtreme(weightfunction,n=length(BL)-1,D=D,mini=FALSE,msg=verbose)
    D0 <- which.max(tmp$extreme-pen)
    lhvalue<-max(tmp$extreme-pen)-n*log(b-a)-n*log(n)
  }
  if (penalty=="aic") {
    pen <- cont2$alpha*(0:(D-1))
    # compute dynamically the max likelihood for irregular histograms with
    # d bins, d varying between 1 and D
    tmp <- DynamicExtreme(weightfunction,n=length(BL)-1,D=D,mini=FALSE,msg=verbose)
    D0 <- which.max(tmp$extreme-pen)
    lhvalue<-max(tmp$extreme-pen)-n*log(b-a)-n*log(n)
  }
  if (penalty=="bic") {
    pen <- cont2$alpha*log(n)*(0:(D-1))
    # compute dynamically the max likelihood for irregular histograms with
    # d bins, d varying between 1 and D
    tmp <- DynamicExtreme(weightfunction,n=length(BL)-1,D=D,mini=FALSE,msg=verbose)
    D0 <- which.max(tmp$extreme-pen)
    lhvalue<-max(tmp$extreme-pen)-n*log(b-a)-n*log(n)
  }
  if (penalty=="cv") {
    # compute dynamically the cv scores for irregular histograms with
    # d bins, d varying between 1 and D
    tmp <- DynamicExtreme(weightfunction,n=length(BL)-1,D=D,mini=TRUE,msg=verbose)
    D0 <- which.min(tmp$extreme)
  }
    
  if (verbose) message(paste("- Number of bins chosen: ",D0,".\n\n",sep=""))
  bounds <- DynamicList(tmp$ancestor,BL,D0)

  # transform back to original range
  y<-a+(b-a)*y
  
  # create histogram
  H <- hist(y,breaks=a+(b-a)*bounds,right=right,plot=FALSE)
  
  # Bugfix: Name of y-var gets lost above - reset it.
  H$xname = yvarname

  # Plot
  if (plot)
    plot( H, freq=FALSE )

  return(list(H=H,lhvalue=lhvalue))
}

