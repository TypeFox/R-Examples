#' Basic summary statistics for ``big.matrix'' objects
#' 
#' @description Functions operate on columns of a 
#' \code{\link[bigmemory]{big.matrix}} object
#' @aliases colmin,min,colmax,max,colprod,prod,colsum,sum,colrange,range,colmean,mean,colvar,var,colsd,sd,colna,summary
#' @param x a \code{\link[bigmemory]{big.matrix}} object.
#' @param object a \code{\link[bigmemory]{big.matrix}} object.
#' @param cols a scalar or vector of column(s) to be summarized.
#' @param na.rm if \code{TRUE}, remove \code{NA} values before summarizing.
#' @param \dots options associated with the correspoding default \R function.
#' @useDynLib biganalytics
#' @details These functions essentially apply summary functions to each 
#' column (or each specified column) of the 
#' \code{\link[bigmemory]{big.matrix}} in turn.
#' @return For \code{colrange}, a matrix with two columns and 
#' \code{length(cols)} rows; column 1 contains the minimum, and column 2 
#' contains the maximum for that column.  The other functions return vectors 
#' of length \code{length(cols)}.
#' @rdname summary.stat
#' @export
#' @examples
#' x <- as.big.matrix(
#'   matrix( sample(1:10, 20, replace=TRUE), 5, 4,
#'           dimnames=list( NULL, c("a", "b", "c", "d")) ) )
#' x[,]
#' mean(x)
#' colmean(x)
#' colmin(x)
#' colmin(x, 1)
#' colmax(x)
#' colmax(x, "b")
#' colsd(x)
#' colrange(x)
#' range(x)
#' colsum(x)
#' colprod(x)
setGeneric('colmin', function(x, cols=NULL, na.rm=FALSE)
  standardGeneric('colmin'))

#' @rdname summary.stat
#' @export
setMethod('colmin', signature(x='big.matrix'),
  function(x, cols=NULL, na.rm=FALSE) {
    thistype = bigmemory:::getCType(x)
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    #if (is.null(cols)) cols = 1:ncol(x)
    #if (is.character(cols)) cols <- mmap(cols, colnames(x))
    ret = .Call("CMinColmain", as.integer(thistype), x@address, 
      as.double(cols), na.rm, PACKAGE="biganalytics")
    if (!is.null(colnames(x))) 
      names(ret) = colnames(x)[cols]
    return(ret)
  })

#' @rdname summary.stat
#' @export
setMethod("min", signature="big.matrix",
  function(x, ..., na.rm=FALSE) {
    return(min(colmin(x, ..., na.rm=na.rm)))
  })

#' @rdname summary.stat
#' @export
setGeneric('colmax', function(x, cols=NULL, na.rm=FALSE)
  standardGeneric('colmax'))

# TODO: Can this be optimized to go through a set of rows only once?
#' @rdname summary.stat
#' @export
setMethod('colmax', signature(x='big.matrix'),
  function(x, cols=NULL, na.rm=FALSE) {
    thistype = bigmemory:::getCType(x)
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    ret = .Call("CMaxColmain", as.integer(thistype), 
      x@address, as.double(cols), na.rm,
      PACKAGE="biganalytics")

    if (!is.null(colnames(x))) 
      names(ret) = colnames(x)[cols]
    return(ret)
  })

#' @rdname summary.stat
#' @export
setMethod("max", signature="big.matrix",
  function(x, ..., na.rm=FALSE)
  {
		return(max(colmax(x, ..., na.rm=na.rm)))
  })

#' @rdname summary.stat
#' @export
setGeneric('colprod', function(x, cols=NULL, na.rm=FALSE)
  standardGeneric('colprod'))

# TODO: Can this be optimized to go through a set of rows only once?
#' @rdname summary.stat
#' @export
setMethod('colprod', signature(x='big.matrix'),
  function(x, cols=NULL, na.rm=FALSE) {
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    #if (is.null(cols)) cols = 1:ncol(x)
    #if (is.character(cols)) cols <- mmap(cols, colnames(x))
    thistype = bigmemory:::getCType(x)
    ret = .Call("CProdColmain", as.integer(thistype), x@address, 
      as.double(cols), na.rm, PACKAGE="biganalytics")
    if (!is.null(colnames(x))) 
      names(ret) = colnames(x)[cols]
    return(ret)
  })

#' @rdname summary.stat
#' @export
setMethod("prod", signature="big.matrix",
  function(x, ..., na.rm=FALSE) {
    return(prod(colprod(x, ..., na.rm=na.rm)))
  })

#' @rdname summary.stat
#' @export
setGeneric('colsum', function(x, cols=NULL, na.rm=FALSE)
  standardGeneric('colsum'))

#' @rdname summary.stat
#' @export
setMethod('colsum', signature(x='big.matrix'),
  function(x, cols=NULL, na.rm=FALSE) {
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    thistype = bigmemory:::getCType(x)
    ret = .Call("CSumColmain", as.integer(thistype), x@address, 
      as.double(cols), na.rm, PACKAGE="biganalytics")
    if (!is.null(colnames(x))) 
      names(ret) = colnames(x)[cols]
    return(ret)
  })

#' @rdname summary.stat
#' @export
setMethod("sum", signature="big.matrix",
  function(x, ..., na.rm=FALSE) {
    return(sum(colsum(x, ..., na.rm=na.rm)))
  })


#' @rdname summary.stat
#' @export
setGeneric('colrange', function(x, cols=NULL, na.rm=FALSE)
  standardGeneric('colrange'))

#' @rdname summary.stat
#' @export
setMethod('colrange', signature(x='big.matrix'),
  function(x, cols=NULL, na.rm=FALSE) {
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    ret = matrix(c(colmin(x,cols=cols,na.rm=na.rm), 
      colmax(x,cols=cols,na.rm=na.rm)), ncol=2)
    colnames(ret) = c('min', 'max')
    if (!is.null(colnames(x))) rownames(ret) = colnames(x)[cols]
    return(ret)
  })

#' @rdname summary.stat
#' @export
setMethod("range", signature="big.matrix",
  function(x, ..., na.rm=FALSE)
  {
    rangeMat = colrange(x, ..., na.rm=na.rm)
    return(c(min(rangeMat[,1]), max(rangeMat[,2])))
  })

#' @rdname summary.stat
#' @export
setGeneric('colmean', function(x, cols=NULL, na.rm=FALSE) 
  standardGeneric('colmean'))

#' @rdname summary.stat
#' @export
setMethod('colmean', signature(x='big.matrix'),
  function(x, cols=NULL, na.rm=FALSE) 
  {
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    #if (is.null(cols)) cols=1:ncol(x)
    #if (is.character(cols)) cols <- mmap(cols, colnames(x))
    thistype = bigmemory:::getCType(x)
    ret = .Call("CMeanColmain", as.integer(thistype), x@address, 
      as.double(cols), na.rm, PACKAGE="biganalytics")
    if (!is.null(colnames(x))) 
      names(ret) = colnames(x)[cols]
    return(ret)
  })

#' @rdname summary.stat
#' @export
setMethod('mean', signature(x="big.matrix"),
  function(x, ...)
  {
    return(mean(colmean(x, ...)))
  })

#' @rdname summary.stat
#' @export
setGeneric('colvar', function(x, cols=NULL, na.rm=FALSE) 
  standardGeneric('colvar'))

#' @rdname summary.stat
#' @export
setMethod('colvar', signature(x='big.matrix'),
  function(x, cols=NULL, na.rm=FALSE)
  {
    if (!is.big.matrix(x)) stop("Unknown type.")
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    #if (is.null(cols)) cols = 1:ncol(x)
    #if (is.character(cols)) cols <- mmap(cols, colnames(x))
    thistype = bigmemory:::getCType(x)
    ret = .Call("CVarColmain", as.integer(thistype), x@address, 
      as.double(cols), na.rm, PACKAGE="biganalytics")
    if (!is.null(colnames(x))) 
      names(ret) = colnames(x)[cols]
    return(ret)
  })

#' @rdname summary.stat
#' @export
setGeneric('colsd', function(x, cols=NULL, na.rm=FALSE)
  standardGeneric('colsd'))

#' @rdname summary.stat
#' @export
setMethod('colsd', signature(x='big.matrix'),
  function(x, cols=NULL, na.rm=FALSE)
  {
    return(sqrt(colvar(x, cols=cols, na.rm=na.rm)))
  })

#' @rdname summary.stat
#' @export
setGeneric('colna', function(x, cols=NULL) standardGeneric('colna'))

#' @rdname summary.stat
#' @export
setMethod('colna', signature(x='big.matrix'),
  function(x, cols=NULL)
  {
    cols <- bigmemory:::cleanupcols(cols, ncol(x), colnames(x))
    #if (is.null(cols)) cols = 1:ncol(x)
    #if (is.character(cols)) cols <- mmap(cols, colnames(x))
    cols = as.double(cols)
    if (max(cols) > ncol(x) | min(cols) < 1) stop("Invalid columns")
    ret = c()
    for (col in cols) ret = c(ret, .Call('ColCountNA', x@address, col, 
                                         PACKAGE="biganalytics"))
    if (!is.null(colnames(x))) names(ret) = colnames(x)[cols]
    return(ret)
  })

#setMethod('summary',
#  signature(object='big.matrix'),
#  function(object)
#  {
#    rows = 1:ncol(object)
#    cn = c('min', 'max', 'mean', "NAs")
#    s = matrix(NA, ncol = length(cn), nrow = length(rows))
#    colnames(s) = cn
#    rownames(s) = colnames(object)
#    s[,'min'] = colmin(object, rows, na.rm=TRUE)
#    s[,'max'] = colmax(object, rows, na.rm=TRUE)
#    s[,'mean'] = colmean(object, rows, na.rm=TRUE)
#    s[,"NAs"] = colna(object, rows)
#    tab=as.table(s)
#    return(tab)
#  })

#' @rdname summary.stat
#' @export
setMethod('summary',
  signature(object='big.matrix'),
  function(object)
  {
    rows <- 1:ncol(object)
    cn <- c('min', 'max', 'mean', 'NAs')
    s <- matrix(NA, ncol = length(cn), nrow = length(rows))
    colnames(s) <- cn
    rownames(s) <- colnames(object)
    for (i in rows) {
      s[i,'min'] <- colmin(object, i, na.rm=TRUE)
      s[i,'max'] <- colmax(object, i, na.rm=TRUE)
      s[i,'mean'] <- colmean(object, i, na.rm=TRUE)
      s[i,'NAs'] <- colna(object, i)
    }
    tab <- as.table(s)
    return(tab)
  })



