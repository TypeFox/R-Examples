###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: dtw.R 388 2015-05-19 19:09:08Z tonig $
#                                                             #
###############################################################


##
## Frontend stuff, including coercing shorthand types
##


`dtw` <-
function(x, y=NULL,
         dist.method="Euclidean",
         step.pattern=symmetric2,
         window.type="none",
         keep.internals=FALSE,
         distance.only=FALSE,
         open.end=FALSE,
         open.begin=FALSE,
         ... ) {

  lm <- NULL;


  
  ## if matrix given
  if(is.null(y)) {
      if(!is.matrix(x)) 
        stop("Single argument requires a global cost matrix");
    
      lm <- x;
  } else if(is.character(dist.method)) {
      ## two timeseries or vectors given
      ## as.matrix coerces ts or mts to matrices
      x <- as.matrix(x);
      y <- as.matrix(y);
      lm <- proxy::dist(x,y,method=dist.method);
  } else if(is.function(dist.method)) {
      stop("Unimplemented");
  } else {
      stop("dist.method should be a character method supported by proxy::dist()");
  }


  

  ## Now we have a function
  wfun<-.canonicalizeWindowFunction(window.type);
  

  ## Now we have a step pattern
  dir<-step.pattern;
  norm <- attr(dir,"norm");


  ## Warn for obsolete constructs
  if(! is.null(list(...)$partial) ) {
    warning("Argument `partial' is obsolete. Use `open.end' instead");
    open.end <- TRUE;
  }



  ## shorthand names
  n <- nrow(lm);
  m <- ncol(lm);

  
  ## For open-begin alignment:
  if (open.begin) {
    
    ##  ensure proper normalization
    if(is.na(norm) || norm != "N") {
      stop("Open-begin requires step patterns with 'N' normalization (e.g. asymmetric, or R-J types (c)). See papers in citation().");
    }

    ## prepend a null row
    lm <- rbind(0,lm);
    np <- n+1;

    ##  pre-initialize elements in the cumulative cost matrix
    precm <- matrix(NA,nrow=np,ncol=m);
    precm[1,] <- 0;

  } else {
    precm <- NULL;
    np <- n;
  }

  
  ## perform the computation
  gcm <- globalCostMatrix(lm, step.matrix=dir,
                          window.function=wfun,
                          seed=precm, ...);


  ## remember size
  gcm$N <- n;
  gcm$M <- m;

  ## remember  call
  gcm$call <- match.call();
  gcm$openEnd <- open.end;
  gcm$openBegin <- open.begin;
  gcm$windowFunction <- wfun;

  ## last row (misnamed), normalized
  lastcol <- gcm$costMatrix[np,];

  if(is.na(norm)) {
      # NO-OP
  } else if(norm == "N+M") {
      lastcol <- lastcol/(n+(1:m));
  } else if(norm == "N") {
      lastcol <- lastcol/n;
  } else if(norm == "M") {
      lastcol <- lastcol/(1:m);
  }

  
  ## for complete alignment
  gcm$jmin <- m;

  
  ## for open-end alignment: normalize
  if (open.end) {
    if(is.na(norm)) {
      stop("Open-end alignments require normalizable step patterns");
    }
    gcm$jmin <- which.min(lastcol);
  }

  ## result: distance 
  gcm$distance <- gcm$costMatrix[np,gcm$jmin];

  ## alignment valid?
  if(is.na(gcm$distance)) {
    stop("No warping path exists that is allowed by costraints"); 
  }
  
  
  ## normalized distance
  if(! is.na(norm)) {
      gcm$normalizedDistance <- lastcol[gcm$jmin];
  } else {
      gcm$normalizedDistance <- NA;
  }

  
  if(!distance.only) {
    ## perform the backtrack
    mapping <- backtrack(gcm);
    gcm <- c(gcm,mapping);    ## Add the properties to gcm
  }


  ## open-begin: discard first elements
  if(open.begin) {
    gcm$index1 <- gcm$index1[-1]-1;
    gcm$index1s <- gcm$index1s[-1]-1;
    gcm$index2 <- gcm$index2[-1];
    gcm$index2s <- gcm$index2s[-1];
    lm <- lm[-1,];
    gcm$costMatrix <- gcm$costMatrix[-1,];
    gcm$directionMatrix <- gcm$directionMatrix[-1,];
  }


  ## don't keep internals: delete sizey intermediate steps 
  if(!keep.internals) {
      gcm$costMatrix<-NULL;
      gcm$directionMatrix<-NULL;
  } else {
  ## keep internals: add data
      gcm$localCostMatrix <- lm;
      if(! is.null(y)) {
          gcm$query <- x;
          gcm$reference <- y;
      }
  }


  ## if a dtw object is to be sponsored:
  class(gcm) <- "dtw";
  return(gcm);
}


##############################
## OO class check
is.dtw <- function(d) {
    return(inherits(d,"dtw"));
}



##############################
## OO print method
print.dtw <- function(x,...) {
  head <- "DTW alignment object\n";
  size <- sprintf("Alignment size (query x reference): %d x %d\n",x$N,x$M);
  call <- sprintf("Call: %s\n",deparse(x$call));
  cat(head,size,call);
}




## Replace  char window.type  with appropriate
## windowing FUNCTION 

.canonicalizeWindowFunction <- function(w) {
  if(is.function(w)) {
    return(w);
  }

  # else 
  wt<-pmatch(w,c("none","sakoechiba","itakura","slantedband"));
  if(is.na(wt)) {
    stop("Ambiguous or unsupported char argument for window.type");
  } 

  wfun<-switch(wt,
	noWindow,
	sakoeChibaWindow,
	itakuraWindow,
	slantedBandWindow);

  return(wfun);
}


