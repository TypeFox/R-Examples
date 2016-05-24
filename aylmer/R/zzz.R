"no.of.boards" <- function(x,n=1e5){
  stopifnot(is.matrix(x))
  if(n>0){
    n <- n+1
  }
  jj <- c("numboards", .Cargs(x), ans=as.integer(n), PACKAGE="aylmer")
  out <- do.call(".C",jj)$ans
  if(out == n){
    warning("Full number of boards exceeds n")
  }
return(out)
}  

"allboards" <- function(x, n=1e5, func=NULL){
  if(!is.null(func)){
    return(apply(Recall(x,n=n),3,func))
  }
  stopifnot(is.matrix(x))
  n <-  no.of.boards(x,n=n)

  flash <- 
    c("allboards",
      .Cargs(x),
      list(out=integer(length(x) * n)),
      as.integer(n),
      PACKAGE="aylmer")
  out <- do.call(".C", flash)$out
  dim(out) <- c(dim(x) , n)
  if(!is.null(dimnames(x))){
    dimnames(out) <- c(dimnames(x), list(NULL))
  }
  jj <- marginals(x)$na
  out[cbind(kronecker(jj , rep(1,n)),rep(seq_len(n) , nrow(jj)))] <- NA
  return(out)
}

"allprobs" <- function(x, n=1e5, normalize=TRUE, give.log=FALSE, use.C=TRUE){
  stopifnot(is.matrix(x))
  n <- no.of.boards(x, n=n)
  if(use.C){
    flash <-
      c("allboardprobs",
        .Cargs(x),
        list(ans=double(n)),
        list(as.integer(n)), 
        PACKAGE="aylmer"
        )
    out <- do.call(.C, flash)$ans
  } else {
    myprob <- function(x){exp(-sum(lfactorial(x),na.rm=TRUE)) }
    jj <- allboards(x, n=n)
    out <- apply(jj,3,myprob)
  }
  if(give.log){
    if(normalize){
      out <- out-max(out)
    }
  } else {
    out <- exp(out)
    if(normalize){
      return(out/sum(out))
    }
  }
  return(out)
}

".facbrob" <- function(x){
  if(is.complex(x)){
    stop("not implemented for complex values")
  }
  test <- x>170
  out <- as.brob(x*0)
  jj <- x[test]+1
  if(any(test)){
    out[test] <-  sqrt(2*pi/jj) * (brob(-1)*( (jj) + 1/(12*jj-0.1/jj)))^jj
  }
  out[!test] <- as.brob(factorial(x[!test]))
  return(out)
}

"prob" <- function(x, give.log=TRUE, use.brob=FALSE){
  stopifnot(is.matrix(x))
  x <- as.vector(x[!is.na(x)])
  if(any(x<0)){
    if(give.log){
      out <- -Inf
    } else {
      out <- 0
    }
      if(use.brob){
        return(as.brob(out))
      } else {
        return(out)
      }
  }
  if(give.log){
    if(use.brob){
      out <- -sum(log(.facbrob(x)))
    } else {
      out <- -sum(lfactorial(x))
    }
  } else {
    if(use.brob){
      out <- 1/prod(.facbrob(x))
    } else {
      out <- 1/prod(factorial(x))
    }
  }
  return(out)

  ## Following code calls C++ routine:

  ## flash <- c(name="prob", .Cargs(x), list(as.integer(x)),
  ##    list(ans=double(1)), PACKAGE="aylmer")
  ##    return(exp(do.call(.C,flash)$ans))

}
