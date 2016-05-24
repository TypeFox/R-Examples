MgSplitMethods <- function(splitMethod,B,N,M,k){
  ReName <- NULL
  k <- as.numeric(substring(grep("^cv[0-9]+$",splitMethod,value=TRUE,ignore.case=TRUE),3))
  if (length(k)==0) k <- NULL
  if (is.null(k)){
    re.noinf <- length(grep("noinf",splitMethod,value=FALSE,ignore.case=TRUE))>0
    re.boot <- length(grep("boot|plain",splitMethod,value=FALSE,ignore.case=TRUE))>0
    re.bootcv <- length(grep("bootcv|out.of.bag|out-of-bag|b0|bootcv",splitMethod,value=FALSE,ignore.case=TRUE))>0
    re.632 <- length(grep("632",splitMethod,value=FALSE,ignore.case=TRUE))>0
    re.plus <- length(grep("plus|\\+",splitMethod,value=FALSE,ignore.case=TRUE))>0
    ## splitMethod <- match.arg(splitMethod,c("none","plain","bootcv","boot632","boot.632","boot632plus","boot.632plus","noinf"))
    if (re.noinf==TRUE){splitMethod <- "noinf"; ReName <- "no information"}
    else if (re.bootcv==TRUE){splitMethod <- "bootcv"; ReName <- "bootcv"}
    else
      if (re.boot==TRUE){
        if (re.632==TRUE){
          if (re.plus==TRUE){splitMethod <- "boot632plus"; ReName <- ".632+"}
          else{splitMethod <- "boot632";ReName <- ".632"}
        }
        else{stop("SplitMethod boot not supported.");splitMethod <- "plain"; ReName <- "bootstrap"}
      }
    if (is.null(ReName)) {splitMethod <- "noSplitMethod"; ReName <- "full data"}
  }
  else{
    splitMethod <- "crossval"
    ReName <- paste(k,"fold cross-validation",sep="-")
  }
  if (missing(M)) M <- N
  stopifnot(M>0 && M<=N) 
  subsampling <- M!=N
  ## if (missing(na.accept)) na.accept <- M/10
  if (splitMethod=="noSplitMethod"|| splitMethod=="noinf") {
    B <- 0
  }
  else{
    if (missing(B)){
      if (length(k)>0) B <- 1 # repeat k-fold CrossVal ones
      else B <- 100  # either `plain' or `Bootcv'
    }
    else if (B==0) stop("No. of crossvals must be a positive integer.")
  }
  if (length(k)>0){
    CrossvalIndex <- do.call("cbind",lapply(1:B,function(b){sample(rep(1:k,length.out=N))}))
  }
  else{
    if (splitMethod %in% c("boot632plus","bootcv","boot632","plain")){
      CrossvalIndex <- matrix(sapply(1:B,function(b){sort(sample(1:N,size=M,replace=!subsampling))}),nrow=M,ncol=B)
      colnames(CrossvalIndex) <- paste("Boot",1:B,sep=".")
    }
    else{
      CrossvalIndex <- NULL
    }
  }
  out <- list(name=ReName,internal.name=splitMethod,index=CrossvalIndex,k=k,B=B,M=M,N=N)
  class(out) <- "MgSplitMethods"
  out
}

