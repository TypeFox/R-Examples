maximalSets <- function(setlist, index=FALSE){
 if (length(setlist)<=1){
    if (index)
      return(1)
    else
      return(setlist)
  }

   lenx     <- c(lapply(setlist,length), recursive=TRUE)
   ooo      <- order(lenx, decreasing=TRUE)
   setlist2 <- setlist[ooo]
   ends     <- cumsum(c( lapply(setlist2,length), recursive=TRUE ))
   iii<-.C("C_maxset",
            setlist=as.character(c(setlist2, recursive=TRUE)),
            ends=ends, nset=length(setlist2), ans=integer(length(setlist2))
           , PACKAGE="gRbase")$ans
   iii <- iii[order(ooo)]

  if (index){
    iii
  } else {
    setlist[iii==1]
  }
}

minimalSets <- function(setlist, index=FALSE){
 if (length(setlist)<=1){
    if (index)
      return(1)
    else
      return(setlist)
  }

   lenx     <- c(lapply(setlist,length), recursive=TRUE)
   ooo      <- order(lenx, decreasing=FALSE)
   setlist2 <- setlist[ooo]
   ends     <- cumsum(c( lapply(setlist2,length), recursive=TRUE ))
   iii<-.C("C_minset",
            setlist=as.character(c(setlist2, recursive=TRUE)),
            ends=ends, nset=length(setlist2), ans=integer(length(setlist2))
           , PACKAGE="gRbase")$ans
   iii <- iii[order(ooo)]

  if (index){
    iii
  } else {
    setlist[iii==1]
  }
}



## A function to remove redundant generators.  If maximal=T, returns
## the maximal generators, if =F, the minimal generators.
## Can be speeded up if the as.character part can be avoided...

removeRedundant <- function(setlist, maximal=TRUE, index=FALSE){
  if (maximal)
    maximalSets(setlist, index)
  else
    minimalSets(setlist, index)
}



## Is x contained in any vector in setlist;
is.insetlist <- function(x, setlist, index=FALSE){
  isin(setlist, x, index)
}

isin <- function(setlist, x, index=FALSE){
  len.setlist <- length(setlist)
  if (len.setlist==0){
    if (index)
      return(0)
    else
      return(FALSE)
  }

  if (len.setlist < 1){
    if (index)
      return(rep(1,length(x)))
    else
      return(TRUE)
  }

  ll    <- cumsum(c( lapply(setlist,length), recursive=TRUE ))
  iii<-.C("C_isin",
          as.character(x), length(x),
          as.character(c(setlist,recursive=TRUE)), ll, len.setlist,
          ans=integer(len.setlist)
          , PACKAGE="gRbase")$ans

  if (index) {
    return(iii)
  } else {
    return(any(iii))
  }
}

## Faster versions of 'standard R functions'
## FIXME: subsetof : Rcpp implementation
is.subsetof <- function(x, set){
  #all(.Internal(match(x,  set,  NA_integer_,  NULL))>0)
  all(match(x,set)>0)
}

subsetof <- function(x, y){
  #all(.Internal(match( x, y, 0, NULL))>0)
  all(match(x,y,0)>0)
}

