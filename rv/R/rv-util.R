
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## UTILITIES in the package rv
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  .show(x)
##  .matrix.ragged.array(x)
##  .expand.as.matrix(M)
##  .listByName(x)
##  .permut
##  .dim.index
##  .bracket
##  .impute.by.name
##  .setDimensionByName
##  .rowmax
##  .slice
##  .nodups
##
## ========================================================================
## .show - show a given variable and its value (useful for debugging purposes)
## ========================================================================

.myc <- function (x) {
  paste("c(",paste(x,collapse=","),")", sep="")
}

.show <- function (x) {
  cat(deparse(substitute(x)),"<-", .myc(x), "\n")
}


## ========================================================================
## .matrix.ragged.array - make a matrix out of a ragged array
## ========================================================================
## not used now

.matrix.ragged.array <- function(x) {
  ### DEBUG: need something more efficient here!
  ncol <- max(sapply(x, length))
  s <- array(NA, c(length(x),ncol))
  for (i in seq(length(x))) {
    y <- x[[i]]
    s[i, seq(length(y))] <- y
  }
  s
}

## ========================================================================
## .expand.as.matrix - build a large matrix out of smaller matrices or scalars, diagonally
## ========================================================================
## Useful when building a large covariance matrix from a series of small ones.
##
## Parameters:  A rv or a list.
## Required:    none
## History:     2004-06-  : 
##

.expand.as.matrix <- function(M) ## NOEXPORT
{
  ## M is a rv or a list.
  dims <- c(0,0)
  for (i in 1:length(M)) {
    dims <- dims + if (is.null(dim(M[[i]]))) rep(length(M[[i]]),2) else dim(M[[i]])
  }
  m <- matrix(0,nrow=dims[1],ncol=dims[2])
  col.offset <- 0
  row.offset <- 0
  for (i in 1:length(M)) {
    x <- if (is.null(dim(M[[i]]))) diag(M[[i]],nrow=length(M[[i]])) else M[[i]]
    d <- dim(x)
    m[row.offset+1:d[1],col.offset+1:d[2]] <- x
    row.offset <- row.offset + d[1]
    col.offset <- col.offset + d[2]
  }
  m
}



## .permut: index of multidimensional array.
## 
## Incidentally we can get the index of the vector by
##  z <- .permut(dims)
##  index <- 1 + (z-1) %*% cumprod(z-1)
##  so index is a vector 1:prod(dims).
##

## ========================================================================
## .permut  -  make  permutations c(1,1) ... c(m,n)
## ========================================================================
## Example: 
##   .permut(c(2,3)) gives a matrix
##        [,1] [,2]
##   [1,]    1    1
##   [2,]    2    1
##   [3,]    1    2
##   [4,]    2    2
##   [5,]    1    3
##   [6,]    2    3
##

.permut <- function(dims) ## NOEXPORT
{
  tl <- prod(dims)
  co <- NULL
  p <- 1
  for (i in 1:length(dims)) {
    d <- dims[i]
    x <- rep(1:d, each=p, length.out=tl)
    co <- cbind(co, x)
    p <- p*d
  }
  dimnames(co) <- NULL
  co
}

## ========================================================================
## .dimindex  -  make a list of indices of a matrix
## ========================================================================
##

.dimindex <- function(x) ## NOEXPORT
{
  if (is.null(dim(x))) {
    ix <- .permut(length(x))
  } else {
    ix <- .permut(dim(x))
  }
  ixt <- paste('[',apply(ix, 1, paste, collapse=','),']',sep='')
  ixt
}

## ========================================================================
## .leftprepend  -  prepend spaces for short lines
## ========================================================================
##

.leftadjust <- function (x) ## NOEXPORT
{
  m.x <- max(nchar(x))
  a <- (nchar(x)<m.x)
  if (any(a)) {
    d <- m.x - nchar(x)
    x[a] <- paste(' ', x[a], sep='') ## only ONE single space
  }
  x
}

## ========================================================================
## .dim.index  -  make a list of indices of a matrix
## ========================================================================
##

.dim.index <- function(x, leftadjust=TRUE) ## NOEXPORT
{
  if (is.null(dim(x)))
    ix <-.permut(length(x))
  else 
    ix <- .permut(dim(x))
  ixt <- paste('[',apply(ix, 1, paste, collapse=','),']',sep='')
  if (leftadjust) .leftadjust(ixt) else ixt
}

## ========================================================================
## .bracket  -  x[[ i ]] but recycle if i > length(x)
## ========================================================================

.bracket <- function(x,k) ## NOEXPORT
{
  l <- length(x)
  if(k<=l) x[[k]] else x[[1+((k-1)%%l)]]
}

## ========================================================================
## .slice  -  apply "[ row ]" to arguments but recycle if row > length(argument)
## ========================================================================
## 
## Description: Like lapply(list(...),"[",row) but recycles.
## Parameters:  
## Required:    none
## History:     2004-06-  : 
##

.slice <- function(..., row) ## NOEXPORT
{
  UseMethod('.slice')
}

.slice.default <- function(..., row) ## NOEXPORT
{
  list(...)
}

.slice.rv <- function(..., row) ## NOEXPORT
{
  v <- sapply(list(...), .bracket, row)
  class(v) <- class(rv())
  return(v)
}

## ========================================================================
## .bracket.indices - extract indices from the brackets
## ========================================================================
## Returns NA's for missing indices

.bracket.indices <- function(x, default.bracket=FALSE) ## already in rv 0.923
{
  if (length(x)<1 || !is.character(x)) return(NULL)
  rgx_bracket <- "^(.*\\[(.*)\\].*)$"
  no.brackets <- (regexpr(rgx_bracket, x)<1)
  y <- sub("^(.*\\[(.*)\\].*)$", "\\2", x)
  y <- gsub(",[[:blank:]]*$", ",NA", y) ## DEBUG: won't work for x[,,]
  if (any(no.brackets)) {
    y[no.brackets] <- if (default.bracket) "1" else "0"
  }
  y <- strsplit(y, "[[:blank:]]*,[[:blank:]]*")
  lapply(y, as.numeric)
}


## ========================================================================
## .indices - return the single-digit indices of components 
## ========================================================================
## x : a list of vectors of indices
## dim. : dimension or length of the target matrix or vector.

.indices <- function (x, dim. = NULL) 
{
    if (is.null(x)) 
        return(numeric(0))
    ld <- length(dim.)
    f <- function(x) {
        lx <- length(x)
        if (lx == 1) 
            return(x)
        if (lx != ld || any(x > dim.)) 
            return(NA)
        1 + sum(c(x - 1, 0) * c(1, cumprod(dim.)))
    }
    pos <- sapply(x, f)
    pos
}



## ========================================================================
## .impute.by.name - impute into a vector using the names attribute
## ========================================================================

.impute.by.name <- function (x, y)
{
  ## impute x <- y using the names of y
  n.y <- names(y)
  name.of.y <- deparse(substitute(y))
  if (length(n.y)<1) stop("the names attribute is required for: ", name.of.y)
  bx <- .bracket.indices(n.y)
  ix <- .indices(bx, dim(x))
  if (length(ix)==1 && ix==0) {
    ## Impute the whole vector.
    return(y)
  }
  na.ix <- is.na(ix)
  if (any(na.ix)) {
    warning("Couldn't figure out index ", n.y[which(na.ix)])
    ix <- ix[!na.ix]
    y  <- y[!na.ix]
    if (length(ix)<1) return(x)
  }
  x[ix] <- y
  if (is.null(dim(x)) || length(dim(x))==1) {
    dim(x) <- NULL
    n.x <- names(x)
    n.x[ix] <- names(y)
    n.x[is.na(n.x)] <- ""
    names(x) <- n.x
  }
  return(x)
}



## ========================================================================
## .setDimensionByName - 
## ========================================================================

.setDimensionByName <- function (x) {
  names.x <- names(x)
  bix <- .bracket.indices(names.x, default.bracket=TRUE)
  if (is.null(bix)) stop("Names of x MUST be set")
  b <- sapply(bix, prod)
  max.ix <- which(b==max(b))[1]
  maxdim <- bix[[max.ix]]
  if (prod(maxdim)<1) stop("Invalid dimension")
  a <- rvarray(NA, maxdim) #### was: array. otherwise .impute.by.name won't work
  new.x <- .impute.by.name(a, x)
  names.new.x <- names(new.x)
  if (length(maxdim)>1) {
    dim(new.x) <- maxdim
  }
  names(new.x) <- names.new.x
  new.x
}



## ========================================================================
## .make.names - make names of the components
## ========================================================================

.make.names <- function(x, name=deparse(substitute(x)))
{
  if (is.null(x) || length(x)<1) return(x)
  name <- make.names(name)
  paste(name, .dim.index(x), sep="")
}

## ========================================================================
## .shortnames - names of the components, with brackets removed
## ========================================================================

.shortnames <- function(x)
{
  na <- names(x)
  if (is.null(na)) return(na)
  unlist(lapply(strsplit(na, "[", fixed=TRUE), "[", 1), use.names=FALSE)
}

## ========================================================================
## make.fullnames - make eventually available!!
## ========================================================================

##.make.fullnames <- function(x)
##{
##  if (is.null(x) || length(x)<1) return(x)
##  name <- make.names(name)
##  paste(name, .dim.index(x), sep="")
##}


.list2array <- function (x, drop=FALSE)
{
  ## name:
  ##   list2array - coerce a list to an array, preserving dimnames
  ## description:
  ##   attempt to coerce the list to an array, preserving dimnames
  ## arguments:
  ##   x : (list) list to be coerced into one array
  ## value:
  ##   a vector if not all dimensions of list elements were equal;
  ##   an array of dimension c(dim(x[[1]]),length(x)),
  ##   if all dimensions of the elements were equal.
  ## author:
  ##   Jouni Kerman
  ## history:
  ##   2006-08-28
  ##
  dm.1 <- length(x)
  if (dm.1<1) return(unlist(x))
  da <- lapply(x, function (x) if (is.null(dim(x))) length(x) else dim(x))
  dm.2 <- da[[1]]
  dn.1 <- names(x)
  dn.2 <- dimnames(x[[1]])
  ## 
  a <- NULL
  for (i in seq(along=x)) {
    a <- c(a, x[[i]])
  }
  if (!all(sapply(da, function (d) all(d==dm.2)))) {
    return(a)
  }
  names(a) <- NULL
  dim(a) <- c(dm.2, dm.1)
  ## If the list component does not have dimnames...
  if (is.null(dn.2)) {
    ## ...and it doesn't have a dimension...
    if (is.null(dim(x[[1]]))) {
      ## ...then use the names for dimensions...
      dn.2 <- names(x[[1]])
    } ## ...or else just give up and use NULL for that dimension names.
    dn.a <- list(dn.2, dn.1)
  } else {
    dn.a <- c(dn.2, list(dn.1))
  }
  if (!is.null(names(dn.2))) {
    name.x <- deparse(substitute(x))
    names(dn.a) <- c(name.x, names(dn.2))
  }
  dimnames(a) <- dn.a
  if (drop) drop(a) else a
}



## .nodups

.nodups <- function (x, keep.latest=TRUE)
{
  ## NAME
  ## .nodups - remove list items with duplicate names (keep latest)
  ## ARGUMENTS
  ##   x : a list
  if (!keep.latest) stop("keep.latest=FALSE Not yet implemented")
  if (length(x)<2) return(x)
  a <- list()
  x <- rev(x) ## use the fact that x[[na]] gets the *first* one
  for (na in names(x)) {
    a[[na]] <- x[[na]] ## multiple assignments may occur, but that's ok.
  }
  return(a)
}


.dimind <- function (x, dim.=NULL, MARGIN=1)
{
  ## dimension indicator (generalization of row, col)
  if (is.null(dim.)) {
    d <- dim(x)
  } else {
    d <- dim.
  }
  if (is.null(d)) {
    return(rep.int(1, length(x)))
  }
  X <- array(NA, dim=d)
  if (length(d)==1) {
    X[] <- 1
  } else if (length(d)==2) {
    X <- if (MARGIN==1) row(X) else if (MARGIN==2) col(X) else stop("subscript out of bounds")
  } else if (MARGIN==1) {
    X[] <- (1:d[MARGIN])
  } else {
    X[] <- rep(1:d[MARGIN], each=prod(d[1:(MARGIN-1)]))
  }
  return(X)
}


.rvmeansd <- function (x, names.=c("mean", "sd", "NAS", "n.sims")) ## for convenience
{
  a <- list(dimnames=dimnames(x), dim=dim(x), names=names(x))
  S <- sims(as.rv(x))
  m <- colMeans(S, na.rm=TRUE)
  ns <- rvnsims(x)
  v <- apply(S, MARGIN=2, var, na.rm=TRUE)
  ####v <- ((colSums(S^2, na.rm=TRUE)-ns*(m^2))/(ns-1))
  v[rvnsims(x)==1] <- 0
  if (any(naS <- is.na(S))) {
    NAS <- (colMeans(naS) * 100)
  } else {
    NAS <- rep.int(0, length(x))
  }
  s <- sqrt(v)## 29.08.2012 v2.2.1: sqrt(na.omit(v))
  attributes(m) <- a
  attributes(s) <- a
  L <- list(mean=m, sd=s, NAS=NAS, n.sims=ns)
  names(L) <- names.
  return(L)
}

## end rv-util.R


