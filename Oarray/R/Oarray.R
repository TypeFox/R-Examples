# Arrays with arbitrary offsets

###### V1.1, 31.05.01
###### V1.3 Minor change to documentation

# Used format() in print.Oarray
# tidied out diagnostic assign in [<-.Oarray, changed default for offset

###### V1.4 Minor change to documentation, after email from Kurt

# Separate internal.Rd file, avoiding use of \synopsis
# Using \details to illustrate use of x[i, j] and x[i, j] <-
# Explicit INDEX file showing main functions
# GPL as licence in DESCRIPTION
# Changed email addresses to Bristol
# Used as.array.default <- base::as.array for making as.array generic [now changed for R >=- 2.8.0 ....rksh]

"Oarray" <-
function(data=NA, dim=length(data), dimnames=NULL, offset=rep(1, length(dim)),
  drop.negative=TRUE)
{
  if(length(offset)==1){offset <- rep(offset,length(dim))}
  
  if (!is.numeric(offset) || length(offset) != length(dim))
    stop("\"offset\" must be numeric vector with same length as \"dim\"")
  if (drop.negative && any(offset < 0))
    stop("Non-negative offsets only")

  robj <- array(data = data, dim = dim, dimnames = dimnames)
  attr(robj, "offset") <- offset
  attr(robj, "drop.negative") <- drop.negative
  class(robj) <- "Oarray"
  robj
}

"as.Oarray" <- function(x, offset=rep(1, length(dim)), drop.negative=TRUE)
{
  x <- as.array(x)
  dim <- dim(x)
  Oarray(x, dim = dim, dimnames = dimnames(x), offset = offset,
    drop.negative = drop.negative)
}


if(getRversion() < "2.8.0"){
  as.array <- function(x,...){UseMethod("as.array")}
  as.array.default <- function(x,...){base::as.array(x,...)}
}

"as.array.Oarray" <- function(x,...)
{
  x <- unclass(x)
  attr(x, "offset") <- NULL
  attr(x, "drop.negative") <- NULL
  NextMethod(x)
}

"is.Oarray" <- function(x)
  inherits(x, "Oarray") && !is.null(attr(x, "offset")) &&
    !is.null(attr(x, "drop.negative"))

# this function takes numeric index sets from the original call
# and maps them using the offset: note that drop=FALSE only works
# if provided as the final argument

".handleTheOffset" <- function(mc, dim, offset, dn)
{
  for (i in seq(along=dim)) {
    ii <- mc[[2+i]]

    if (missing(ii)) next

    if (is.symbol(ii) || is.call(ii))
      ii <- eval.parent(ii, 3)

    if (is.numeric(ii)) {

      if (!dn || all(ii>=0))
        ii <- ifelse(ii>=offset[i], ii - offset[i] + 1, dim[i]+1)
      else {
        if (all(ii <= -offset[i]))
          ii <- ii + offset[i] - 1
        else stop("subscript out of bounds")
      }

      mc[[2+i]] <- ii
    }
  }
  mc
}

"[.Oarray" <- function(x, ...)
{
  mc <- match.call()
  k <- length(mc)
  offset <- attr(x, "offset")
  dn <- attr(x, "drop.negative")
  dim <- dim(x)

  if( k==3 ){
    if(mc[[3]] == ""){
      return(as.array(x))
    } 
    args <- list(...)
    index <- args[[1]]
    if(is.logical(index)){
      return(as.array(x)[index])
    }
    if(is.matrix(index)){
      return(as.array(x)[1+sweep(index,2,offset)])
    }
  }

  if (k < 2+length(dim))
    stop("incorrect number of dimensions")

  mc <- .handleTheOffset(mc, dim, offset, dn)
  mc[[1]] <- as.name("[")
  mc[[2]] <- as.name("x")
  x <- as.array(x)
  eval(mc, envir=NULL)
}

"[<-.Oarray" <- function(x, ..., value)
{
  mc <- match.call()
  k <- length(mc)
  offset <- attr(x, "offset")
  dn <- attr(x, "drop.negative")
  dim <- dim(x)

  if (k==4){
    if (mc[[3]] == ""){
      return(Oarray(value, dim, dimnames(x), offset, dn))
    }
    args <- list(...)
    index <- args[[1]]
    att <- attributes(x)
    if(is.logical(index)){
      x <- as.array(x)
      x[index] <- value
      attributes(x) <- att
      return(x)
    }
    if(is.matrix(index)){
      x <- as.array(x)
      x[1+sweep(index,2,offset)] <- value
      attributes(x) <- att
      return(x)
    }
  }
  
  if (k < 3+length(dim))
    stop("incorrect number of dimensions")

  mc <- .handleTheOffset(mc, dim, offset, dn)
  mc[[1]] <- as.name("[<-")
  mc[[2]] <- as.name("x")
  x <- as.array(x)
  robj <- eval(mc, envir=NULL)
  Oarray(robj, dim, dimnames(x), offset, dn)
}

"print.Oarray" <-
function(x, ...)
{
  d <- dim(x)
  dn <- dimnames(x)
  if (is.null(dn))
    dn <- vector("list", length(d))
  offset <- attr(x, "offset")
  x <- as.array(x)

  for (i in seq(along=dn))
    if (is.null(dn[[i]])) {
      dn[[i]] <- 0:(d[i]-1) + offset[i]
      if (i==1 || i==2) {
        dn[[i]] <- format(dn[[i]])
        if (i==1)
          dn[[i]] <- paste("[", dn[[i]], ",]", sep="")
        else
          dn[[i]] <- paste("[,", dn[[i]], "]", sep="")
      }
    }
  dimnames(x) <- dn
  NextMethod("print")
}


setOldClass("Oarray")

setMethod("slice.index","Oarray",function(x,MARGIN){
  o <- attr(x,"offset")
  attr <- attributes(x)
  out <- slice.index(as.array(x),MARGIN)+o[MARGIN]-1L
  attributes(out) <- attr
  return(out)
} )


