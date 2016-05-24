check.input<-function(x, ncat, nitem=NULL, nq=48, mxiter=200, iprint=FALSE) {
    
  if (is.null(nitem)){type="raw"} else {type="freqtable"}

  x<-as.matrix(x)
  
  ## input checking & data prep
  if(!ncat%%1==0){
    warning("Non-integers found in input data. These have been rounded to the nearest integer")
    ncat<-round(ncat)
  }
  
  if (type=="raw"){
    myData<-stacked.table(x)
    nitem<-myData$nitem
    nrec<-myData$nrec
    myX<-myData$x
  } else if (type=="freqtable") {
    if(nitem != ncol(x)-1){
      warning("nitem does not match number of columns of input matrix. Number of items set as ncol(x)-1.")      
    }
    nitem<-ncol(x)-1
    nrec<-nrow(x)
    myX<-x
  } else {
    stop("Unknown input type", type)
  }

  if(!all(myX%%1==0) || !nitem%%1==0 || !nrec%%1==0){
    warning("Non-integers found in input data. These have been rounded to the nearest integer")
    myX<-round(myX)
    nitem<-round(nitem)
    nrec<-round(nrec)
  }

  if(min(myX[1:(nitem*nrec)])<0 || max(myX[1:(nitem*nrec)])!=(ncat-1)){
    stop("Categorical items must have a min value of 0 and count up at 1 integer intervals, with a maximum equal to (ncat-1)")
  }

  ## estimation options
    
  if(iprint){
    iprint<-1      
  } else if(!iprint) {
    iprint<-0
  } else {
    warning("iprint not logical, setting to FALSE")
    iprint<-0
  }
    
  nq<-round(nq)
  if(nq<10){
    warning("Small number of quadrature points selected")      
  }
  if(nq<2){
    nq<-3
  }
    
  mxiter<-round(mxiter)
  if(mxiter<0){
    warning("Max iterations selected to be 0 or negative")
    mxiter<-0
  }

  list(nitem=nitem, ncat=ncat, nrec=nrec, myX=myX, N=sum(myX[,ncol(myX)]), nq=nq, mxiter=mxiter, iprint=iprint)

}


check.alphas<-function(alphas, nitem, ncat){
  myOk<-TRUE
  if (is.null(alphas)){
    myOk<-FALSE
  } else if (length(alphas)!=(nitem*(ncat-1))){
    warning("Length of alphas incompatible with nitem and ncat")
    myOk<-FALSE
  } else if(!is.double(alphas)){
    ## TO DO: think about input type for alphas more carefully and whether
    ## this is necessary / does proper check
    warning("Input for alphas not of type double")
    myOk<-FALSE
  }
  myOk
}

check.betas<-function(bVec, nitem) {
  myOk<-TRUE
  if (is.null(bVec)) {
    myOk<-FALSE
  }
  else if (length(bVec)!=nitem) {
    myOk<-FALSE
  }
  myOk
}

check.bounds<-function(startvals, bounds){
    # vector of length 2
    if(!(length(bounds)==2)){
      stop("Input for lower and upper bounds on parameter estimates not a vector of length 2", bounds)
    }
    #lower less than greater
    if(bounds[1]>=bounds[2]){
      stop("Input for upper bound on parameter estimates equal to or less than lower bound", bounds)
    }
    
    # start values w/in bounds
    for (i in startvals){
      if(i<bounds[1]||i>bounds[2]){
        warning("Starting value outside of boundaries for parameter estimates", i, bounds)
      }
    }
    return(bounds)
}

### 2011-05-23
### Utility function that creates contingency table
### from a data matrix and stacks it in a format
### usable by plnmodel
#stacked.table <- function (x) {
#  newX<-as.data.frame(xtabs(~., data=as.data.frame(x)))
#  newX<-subset(newX, subset=newX[,length(newX)]>0)
#  nrec<-nrow(newX)
#  newX<-as.numeric(as.matrix(newX))
#  nitem<-ncol(x)
#  list(x=newX, nitem=nitem, nrec=nrec)
#}

## more efficient version based on code from the ltm package
stacked.table <- function (x) {
  ##cl <- match.call()

  X <- data.matrix(x)
  ## save colnames, get rid of dimnames

  dimnames(X) <- NULL

  ## estimate number of categories?
  ncatg <- apply(X, 2, function (x) if (any(is.na(x))) length(unique(x)) - 1 else length(unique(x)))

  ## number of rows, number of columns
  n <- nrow(X)
  p <- ncol(X)

  ## converts each row to single strings
  pats <- apply(X, 1, paste, collapse = "/")

  ## frequencies of each pattern
  freqs <- table(pats)

  ## number of patterns
  nfreqs <- length(freqs)

  ## frequencies of each pattern as a vector
  obs <- as.vector(freqs)

  ## patterns as one single vector
  X <- unlist(strsplit(cbind(names(freqs)), "/"))
  X[X == "NA"] <- as.character(NA) ## handle missing values
  
  ## pattern section of contingency table
  X <- matrix(as.numeric(X), nfreqs, p, TRUE)

  list(x=cbind(X,obs), nitem=p, nrec=nfreqs)
}
