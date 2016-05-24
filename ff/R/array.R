# Array utilities for ff
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-09-03
# Last changed: 2007-10-25

# source("D:/mwp/eanalysis/ff/R/array.R")

# xx TODO arrprint, print.arrprint (generalize matprint to arrays)

# an ff dimorder (storage layout) 2,3,1 means that data is stored fastest rotating along columns, then layers, slowest is rows
# Implementation: the R side submits dim[dimorder] to the C-Side (dimorder is transparent for the C-Side: C does not know about the dimorder)
# NOTE that therefor you MUST NEVER access the dim-attribute install("Dim") from C !!
# if we submit ff_array[]<-vector, vector is interpreted in dimorder 1,2,3 (by layer, by column, by row, i.e. row rotates fastest) whatever the ff dimorder
# we can change the intepretation of the vector by ff_array[,dimorder=c(2,3,1)]<-vector, now the vector elements are put fastest col then layer then row
# if we submit ff_array[1:length(vector)]<-vector it is simply copied and ff dimorder DOES matter



# not exported, only used for display in vecprint, matprint
legalizeFactor <- function(x){
  if (is.factor(x)){
    a <- attributes(x)
    attributes(x) <- NULL
    if (any(!is.na(x) & x==0L)){
      x <- x + 1L
      a$levels <- c("<0>", a$levels)
    }
    setattributes(x, a) #attributes(x) <- a
  }
  x
}





#! \name{vector2array}
#! \alias{vector2array}
#! \title{ Array: make array from vector }
#! \description{
#!   makes array from vector respecting \option{dim} and \option{dimorder}
#! }
#! \usage{
#! vector2array(x, dim, dimorder = NULL)
#! }
#! \arguments{
#!   \item{x}{ an input vector, recyled if needed }
#!   \item{dim}{ \code{\link{dim}} }
#!   \item{dimorder}{ \code{\link{dimorder}} }
#! }
#! \details{
#! FILLS vector into array of dim where fastest rotating is dim[dimorder[1]], next is dim[dimorder[2]] and so forth.
#! This is a generalization of converting vector to matrix(, byrow=TRUE).
#! NOTE that the result is a ram array always stored in STANDARD dimorder !!!
#! In this usage we sometimes term the dimorder 'bydim' because it does not change the physical layout of the result,
#! rather bydim refers to the dimorder in which to interpret the vector (not the result).
#! In \command{ff}, \command{update} and \command{clone} we have 'bydim' to contrast it from 'dimorder', the latter describing the layout of the file.
#! }
#! \value{
#!   a suitable \code{\link{array}}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{array2vector}}, \code{\link{vectorIndex2arrayIndex}} }
#! \examples{
#!   vector2array(1:12, dim=c(3, 4))               # matrix(1:12, 3, 4)
#!   vector2array(1:12, dim=c(3, 4), dimorder=2:1) # matrix(1:12, 3, 4, byrow=TRUE)
#! }
#! \keyword{ array }
#! \keyword{ data }

vector2array <- function(x, dim, dimorder=NULL){
  if (is.null(dimorder) || identical(dimorder, 1:length(dim))){
    array(x, dim=dim)
  }else{
    aperm(array(x, dim=dim[dimorder]), match(1:length(dim), dimorder))
  }
}


#! \name{array2vector}
#! \alias{array2vector}
#! \title{ Array: make vector from array }
#! \description{
#!   Makes a vector from an array respecting \option{dim} and \option{dimorder}
#! }
#! \usage{
#! array2vector(x, dim = NULL, dimorder = NULL)
#! }
#! \arguments{
#!   \item{x}{ an \code{\link{array}} }
#!   \item{dim}{ \code{\link{dim}} }
#!   \item{dimorder}{ \code{\link{dimorder}} }
#! }
#! \details{
#!  This is the inverse function of \code{\link{vector2array}}.
#!  It extracts the vector from the array by first moving through the fastest rotating dimension dim[dimorder[1]], then dim[dimorder[2]], and so forth
#! }
#! \value{
#!   a vector
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{vector2array}}, \code{\link{arrayIndex2vectorIndex}} }
#! \examples{
#!   array2vector(matrix(1:12, 3, 4))
#!   array2vector(matrix(1:12, 3, 4, byrow=TRUE), dimorder=2:1)
#! }
#! \keyword{ array }
#! \keyword{ data }

array2vector <- function(x, dim=NULL, dimorder=NULL){
  if (is.array(x)){
    if (is.null(dim))
      dim <- dim(x)
    else if (!identical(dim(x), dim))
      stop("dim(x) differs from dim")
  }else{
    if (is.null(dim))
      stop("need is.array(x) or dim")
  }
  if (is.null(dimorder))
    dimorder <- dimorder(x)
  if (is.null(dimorder)|| identical(dimorder, 1:length(dim)))
    as.vector(x)
  else{
    if (is.array(x))
      as.vector(aperm(x, dimorder))
    else
      as.vector(aperm(array(x, dim=dim), dimorder))
  }
}


if (FALSE){
  # regression tests
  stopifnot(identical( vector2array(1:12, 3:4), matrix(1:12, 3, 4) ))
  stopifnot(identical( vector2array(1:12, 3:4, 2:1), matrix(1:12, 3, 4, byrow=TRUE) ))
  stopifnot(identical( vector2array(1:24, 4:2, c(2,1,3)), structure(c(1L, 4L, 7L, 10L, 2L, 5L, 8L, 11L, 3L, 6L, 9L, 12L, 13L, 16L, 19L, 22L, 14L, 17L, 20L, 23L, 15L, 18L, 21L, 24L), .Dim = c(4L, 3L, 2L)) ))

  stopifnot(identical( array2vector(vector2array(1:12, 3:4)), 1:12 ))
  stopifnot(identical( array2vector(vector2array(1:12, 3:4, 2:1), dimorder=2:1), 1:12 ))
  stopifnot(identical( array2vector(vector2array(1:24, 4:2, c(2,1,3)), dimorder=c(2,1,3)), 1:24 ))
}


#! \name{arrayIndex2vectorIndex}
#! \alias{arrayIndex2vectorIndex}
#! \title{ Array: make vector positions from array index }
#! \description{
#!   Make vector positions from a (non-symmetric) array index respecting \option{dim} and \option{dimorder}
#! }
#! \usage{
#! arrayIndex2vectorIndex(x, dim = NULL, dimorder = NULL, vw = NULL)
#! }
#! \arguments{
#!   \item{x}{ an n by m matrix with n m-dimensional array indices }
#!   \item{dim}{ NULL or \code{\link{dim}} }
#!   \item{dimorder}{ NULL or \code{\link{dimorder}} }
#!   \item{vw}{ NULL or integer vector[3] or integer matrix[3,m], see details }
#! }
#! \details{
#!   The fastest rotating dimension is dim[dimorder[1]], then dim[dimorder[2]], and so forth. \cr
#!   The parameters 'x' and 'dim' may refer to a subarray of a larger array, in this case, the array indices 'x' are interpreted as 'vw[1,] + x' within the larger array 'as.integer(colSums(vw))'.
#! }
#! \value{
#!   a vector of indices in \code{1:prod(dim)} (or  \code{1:prod(colSums(vw))})
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{array2vector}}, \code{\link{vectorIndex2arrayIndex}} }
#! \examples{
#!   x <- matrix(1:12, 3, 4)
#!   x
#!   arrayIndex2vectorIndex(cbind(as.vector(row(x)), as.vector(col(x)))
#!   , dim=dim(x))
#!   arrayIndex2vectorIndex(cbind(as.vector(row(x)), as.vector(col(x)))
#!   , dim=dim(x), dimorder=2:1)
#!   matrix(1:30, 5, 6)
#!   arrayIndex2vectorIndex(cbind(as.vector(row(x)), as.vector(col(x)))
#!   , vw=rbind(c(0,1), c(3,4), c(2,1)))
#!   arrayIndex2vectorIndex(cbind(as.vector(row(x)), as.vector(col(x)))
#!   , vw=rbind(c(0,1), c(3,4), c(2,1)), dimorder=2:1)
#! }
#! \keyword{ array }
#! \keyword{ data }

arrayIndex2vectorIndex <- function(x, dim = NULL, dimorder=NULL, vw=NULL){
  if (is.null(vw)){
    if (is.null(dim))
      stop("need 'dim' or 'vw'")
    n <- length(dim)
  }else{
    n <- ncol(vw)
    if (!identical(dim(vw), c(3L, n)))
      stop("dim(vw) must be c(3L, length(dim))")
  }
  if(!is.matrix(x) || ncol(x)!=n)
    stop("array indexing requires a matrix with ncol(matrix)=length(dim(array))")
  if (is.null(dimorder))
    dimorder <- 1:n
  if (is.null(vw)){
    cdim <- cumprod(c(1, dim[dimorder][-n]))
    as.integer(colSums(t(x[,dimorder,drop=FALSE]-1) * cdim) + 1L)
  }else{
    cdim <- cumprod(c(1, (colSums(vw))[dimorder][-n]))
    as.integer(colSums((vw[1,dimorder] + t(x[,dimorder,drop=FALSE]-1)) * cdim) + 1L)
  }
}

#! \name{vectorIndex2arrayIndex}
#! \alias{vectorIndex2arrayIndex}
#! \title{ Array: make array from index vector positions }
#! \description{
#!   make array from index vector positions respecting \option{dim} and \option{dimorder}
#! }
#! \usage{
#! vectorIndex2arrayIndex(x, dim = NULL, dimorder = NULL, vw = NULL)
#! }
#! \arguments{
#!   \item{x}{ a vector of indices in \code{1:prod(dim)} }
#!   \item{dim}{ NULL or \code{\link{dim}} }
#!   \item{dimorder}{ NULL or \code{\link{dimorder}} }
#!   \item{vw}{ NULL or integer matrix[2,m], see details }
#! }
#! \details{
#!   The fastest rotating dimension is dim[dimorder[1]], then dim[dimorder[2]], and so forth. \cr
#!   The parameters 'x' and 'dim' may refer to a subarray of a larger array, in this case, the array indices 'x' are interpreted as 'vw[1,] + x' within the larger array 'vw[1,] + x + vw[2,]'.
#! }
#! \value{
#!   an n by m matrix with n m-dimensional array indices
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{vector2array}}, \code{\link{arrayIndex2vectorIndex}} , \code{\link{symmIndex2vectorIndex}} }
#! \examples{
#!   matrix(1:12, 3, 4)
#!   vectorIndex2arrayIndex(1:12, dim=3:4)
#!   vectorIndex2arrayIndex(1:12, dim=3:4, dimorder=2:1)
#!   matrix(1:30, 5, 6)
#!   vectorIndex2arrayIndex(c(6L, 7L, 8L, 11L, 12L, 13L, 16L, 17L, 18L, 21L, 22L, 23L)
#! , vw=rbind(c(0,1), c(3,4), c(2,1)))
#!   vectorIndex2arrayIndex(c(2L, 8L, 14L, 3L, 9L, 15L, 4L, 10L, 16L, 5L, 11L, 17L)
#! , vw=rbind(c(0,1), c(3,4), c(2,1)), dimorder=2:1)
#!
#!   \dontshow{
#!     # incomplete regression tests
#!     x <- matrix(1:12, 3, 4)
#!     i <- cbind(as.vector(row(x)), as.vector(col(x)))
#!     stopifnot(identical( array(arrayIndex2vectorIndex(i, dim=dim(x)), dim=dim(x)), x ))
#!     stopifnot(identical( vectorIndex2arrayIndex(arrayIndex2vectorIndex(i, dim=dim(x), dimorder=1:2), dim=dim(x)), i ))
#!
#!     y <- vector2array(1:12, c(3,4), 2:1)
#!     i <- cbind(as.vector(row(y)), as.vector(col(y)))
#!     stopifnot(identical( array(arrayIndex2vectorIndex(i, dim=dim(y), dimorder=2:1), dim=dim(y)), y ))
#!     stopifnot(identical( vectorIndex2arrayIndex(arrayIndex2vectorIndex(i, dim=dim(y), dimorder=dimorder(y)), dim=dim(y), dimorder=dimorder(y)), i ))
#!
#!     z <- vector2array(1:24, dim=4:2, dimorder=c(2,1,3))
#!     stopifnot(identical( arrayIndex2vectorIndex(vectorIndex2arrayIndex(z, dim=dim(z), dimorder=c(2,1,3)), dim=dim(z), dimorder=c(2,1,3)), as.vector(z) ))
#!   }
#! }
#! \keyword{ array }
#! \keyword{ data }


vectorIndex2arrayIndex <- function(x, dim=NULL, dimorder=NULL, vw=NULL){
  if (is.null(vw)){
    if (is.null(dim))
      stop("need 'dim' or 'vw'")
    n <- length(dim)
    if (is.null(dimorder))
      dimorder <- 1:n
    cdim <- as.integer(cumprod(dim[dimorder]))
    vw <- matrix(0L, 3L, n)
  }else{
    n <- ncol(vw)
    if (!identical(dim(vw), c(3L, n)))
      stop("dim(vw) must be c(3L, length(dim))")
    if (is.null(dimorder))
      dimorder <- 1:n
    cdim <- as.integer(cumprod(colSums(vw)[dimorder]))
  }

  N <- cdim[n]
  cdim <- c(1L,cdim[-n])

  if(!is.vector(x))
    x <- as.vector(x)

  if (n==1)
    return(cbind(x-vw[1,]))

  if (!is.integer(x))
    x <- as.integer(x - 1L)
  else
    x <- x - 1L

  ret <- NULL
  for (i in n:2){
    ret <- rbind(x%/%cdim[i], ret)
    x <- x%%cdim[i]
  }
  ret <- rbind(x, ret)
  setattr(ret, "dimnames", NULL) #dimnames(ret) <- NULL
  t(ret[match(1:n, dimorder),,drop=FALSE] - vw[1,] + 1L)
}



#! \name{symmIndex2vectorIndex}
#! \alias{symmIndex2vectorIndex}
#! \title{ Array: make vector positions from symmetric array index }
#! \description{
#!   make vector positions from (non-symmetric) array index respecting \option{dim} and \option{fixdiag}
#! }
#! \usage{
#! symmIndex2vectorIndex(x, dim, fixdiag = NULL)
#! }
#! \arguments{
#!   \item{x}{ a matrix[,1:2] with matrix subscripts }
#!   \item{dim}{ the dimensions of the symmetric matrix }
#!   \item{fixdiag}{ NULL assumes free diagonal, any value assumes fixed diagonal }
#! }
#! \details{
#!   With \option{fixdiag = NULL}
#! }
#! \value{
#!   a vector of indices in \code{1:prod(dim(x))}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{arrayIndex2vectorIndex}} }
#! \examples{
#!   symmIndex2vectorIndex(rbind(
#!    c(1,1)
#!   ,c(1,10)
#!   ,c(10,1)
#!   ,c(10,10)
#!   ), dim=c(10,10))
#!   symmIndex2vectorIndex(rbind(
#!    c(1,1)
#!   ,c(1,10)
#!   ,c(10,1)
#!   ,c(10,10)
#!   ), dim=c(10,10), fixdiag=1)
#! }
#! \keyword{ array }
#! \keyword{ data }



# mapping row/col subscripts to position in symmetric matrix
# in R-counting 1..n
symmIndex2vectorIndex <- function(x, dim, fixdiag=NULL # setting to anything but NULL deactivates the diagonal (e.g. 0 as in dist)
){
  one <- as.integer(1)
  two <- as.integer(2)
  r <- x[,1] - one
  c <- x[,2] - one
  n <- dim[[1]]
  if (is.null(fixdiag)){
    i <- as.integer(ifelse(r>c
    , r + c*n - c*(c+one)/two + one
    , c + r*n - r*(r+one)/two + one
    ))
  }else{
    i <- as.integer(ifelse(r>c
    , r + c*(n-one) - c*(c+one)/two
    , c + r*(n-one) - r*(r+one)/two
    ))
    i[r==c] <- NA
  }
  i
}



#! \name{dummy.dimnames}
#! \alias{dummy.dimnames}
#! \title{ Array: make dimnames }
#! \description{
#!   makes standard dimnames from letters and integers (for testing)
#! }
#! \usage{
#! dummy.dimnames(x)
#! }
#! \arguments{
#!   \item{x}{ an \code{\link{array}} }
#! }
#! \value{
#!   a list with character vectors suitable to be assigned as dimnames to x
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{dimnames}} }
#! \examples{
#!   dummy.dimnames(matrix(1:12, 3, 4))
#! }
#! \keyword{ array }
#! \keyword{ data }

dummy.dimnames <- function(x){
  d <- dim(x)
  if (is.null(d))
    return(NULL)
  lapply(1:length(d), function(n)paste(c(letters,LETTERS)[n], 1:d[n], sep=""))
}


#! \name{matcomb}
#! \alias{matcomb}
#! \title{ Array: make matrix indices from row and columns positions }
#! \description{
#!   create matrix indices from row and columns positions
#! }
#! \usage{
#! matcomb(r, c)
#! }
#! \arguments{
#!   \item{r}{ integer vector of row positions }
#!   \item{c}{ integer vector of column positions  }
#! }
#! \details{
#!   rows rotate faster than columns
#! }
#! \value{
#!   a k by 2 matrix of matrix indices where \code{ k = length(r) * length(c) }
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{row}}, \code{\link{col}} , \code{\link{expand.grid}} }
#! \examples{
#!   matcomb(1:3, 1:4)
#!   matcomb(2:3, 2:4)
#! }
#! \keyword{ array }
#! \keyword{ data }


# create all combinations between r and c
# assuming integer input
matcomb <- function(r,c){
  rc <- rep.int(r, length(c))
  cr <- rep.int(c, rep(length(r), length(c)))
  cbind(r=rc,c=cr)
}


#! \name{matprint}
#! \alias{matprint}
#! \alias{print.matprint}
#! \title{ Print beginning and end of big matrix }
#! \description{
#!   Print beginning and end of big matrix
#! }
#! \usage{
#! matprint(x, maxdim = c(16, 16), digits = getOption("digits"))
#! \method{print}{matprint}(x, quote = FALSE, right = TRUE, \dots)
#! }
#! \arguments{
#!   \item{x}{ a \code{\link{matrix}} }
#!   \item{maxdim}{ max number of rows and columns for printing }
#!   \item{digits}{ see \code{\link{format}} }
#!   \item{quote}{ see \code{\link{print}} }
#!   \item{right}{ see \code{\link{print}} }
#!   \item{\dots}{ see \code{\link{print}} }
#! }
#! \value{
#!   a list of class 'matprint' with components
#!   \item{ subscript }{ a list with four vectors of subscripts: row begin, column begin, row end, column end  }
#!   \item{ example }{ the extracted example matrix as.characer including seperators }
#!   \item{ rsep }{ logical scalar indicating whether row seperator is included }
#!   \item{ csep }{ logical scalar indicating whether column seperator is included }
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{vecprint}} }
#! \examples{
#!   matprint(matrix(1:(300*400), 300, 400))
#! }
#! \keyword{ array }
#! \keyword{ print }


matprint <- function(x, maxdim=c(16,16), digits=getOption("digits")){

  d <- dim(x)
  if (any(d==0)){
    return("[empty matrix]")
  }

  maxdim <- pmin(d,maxdim)
  d2 <- maxdim %/% 2
  d1 <- maxdim - d2
  d2 <- ifelse(d1<maxdim, d - d2 + 1, 0)
  rsep <- maxdim[1] < d[1]
  csep <- maxdim[2] < d[2]

  i <- matrix(list(1:d1[1], 1:d1[2], if (d2[1]) d2[1]:d[1], if (d2[2]) d2[2]:d[2]), 2, 2)
  i1 <- c(i[[1,1]],i[[1,2]])
  i2 <- c(i[[2,1]],i[[2,2]])
  m <- x[ i1 , i2, drop=FALSE ]
  if (is.data.frame(m)){
    for (j in seq(along=m))
      m[[j]] <- legalizeFactor(m[[j]])
    m <- as.matrix(m)
  }else{
    m <- legalizeFactor(m)
  }

  # xx circumvent a bug in format of raw
  if (is.raw(m)){
    a <- attributes(m)
    m <- as.character(m)
    attributes(m) <- a
  }else{
    m <- format(m, digits=digits)
  }
  if (is.null(rownames(m)))
    rownames(m) <- paste("[", c(i[[1,1]],i[[1,2]]),",]", sep="")
  if (is.null(colnames(m)))
    colnames(m) <- paste("[,", c(i[[2,1]],i[[2,2]]),"]", sep="")
  if (csep){
    r1 <- cbind( m[       1 :    d1[1],1:d1[2], drop=FALSE]  , if(d2[2])cbind(" "=":", m[       1 :    d1[1],(d1[2]+1):maxdim[2], drop=FALSE] ))
    r2 <- cbind( if(d2[1])m[(d1[1]+1):maxdim[1],1:d1[2], drop=FALSE]  , if(d2[1]&&d2[2])cbind(" "=":", m[(d1[1]+1):maxdim[1],(d1[2]+1):maxdim[2], drop=FALSE] ))
  }else{
    r1 <- cbind( m[       1 :    d1[1],1:d1[2], drop=FALSE]  ,                if(d2[2])m[       1 :    d1[1],(d1[2]+1):maxdim[2], drop=FALSE] )
    r2 <- cbind( if(d2[1])m[(d1[1]+1):maxdim[1],1:d1[2], drop=FALSE]  ,       if(d2[1]&&d2[2])         m[(d1[1]+1):maxdim[1],(d1[2]+1):maxdim[2], drop=FALSE] )
  }
  if (rsep)
    m <- rbind(r1,":"=":",r2)
  else
    m <- rbind(r1,r2)
  ret <- list(subscript=i, example=m, rsep=rsep, csep=csep)
  class(ret) <- "matprint"
  ret
}
print.matprint <- function(x, quote=FALSE, right=TRUE, ...){
  print(x$example, quote=quote, right=right, ...)
}


#! \name{vecprint}
#! \alias{vecprint}
#! \alias{print.vecprint}
#! \title{ Print beginning and end of big vector }
#! \description{
#!   Print beginning and end of big vector
#! }
#! \usage{
#! vecprint(x, maxlength = 16, digits = getOption("digits"))
#!  \method{print}{vecprint}(x, quote = FALSE, \dots)
#! }
#! \arguments{
#!   \item{x}{ a vector }
#!   \item{maxlength}{ max number of elements for printing }
#!   \item{digits}{ see \code{\link{format}} }
#!   \item{quote}{ see \code{\link{print}} }
#!   \item{\dots}{ see \code{\link{print}} }
#! }
#! \value{
#!   a list of class 'vecprint' with components
#!   \item{ subscript }{ a list with two vectors of subscripts: vector begin and vector end }
#!   \item{ example }{ the extracted example vector as.character including seperator }
#!   \item{ sep }{ the row seperator ":" }
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{matprint}} }
#! \examples{
#!   vecprint(10000:1)
#! }
#! \keyword{ print }


vecprint <- function(x, maxlength=16, digits=getOption("digits")){

  d <- length(x)
  maxlength <- min(maxlength, d)
  d2 <- maxlength%/%2
  d1 <- maxlength - d2
  d2 <- ifelse(d1<maxlength, d - d2 + 1, 0L)

  sep <- maxlength[1]<d[1]

  i <- list(1:d1, if (d2) d2:d)
  m <- format(legalizeFactor(x[ unlist(i) ]), digits=digits)
  if (is.null(names(m)))
    names(m) <- paste("[", unlist(i),"]", sep="")

  if (sep){
    m <- c(m[1:d1], if (d2) c(":", m[(d1+1):maxlength]))
  }
  ret <- list(subscript=i, example=m, sep=sep)
  class(ret) <- "vecprint"
  ret
}
print.vecprint <- function(x, quote=FALSE, ...){
  print(x$example, quote=quote, ...)
}


