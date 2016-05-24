"as.rectangular" <- 
function(x)
{
	if(is.rectangular(x))
		x
	else as.data.frame(x)
}

"as.char.rect" <- 
function(x)
{
  x <- as.rectangular(x)
  dn <- c(NROW(x), NCOL(x))
  ret <- if(is.data.frame(x)) sapply(x, function(col)
                                     if(is.factor(col)) as.matrix(col)
                                     else as.character(col)) else 
  as.character(x)
  if(dn[2] != NCOL(ret))
    dim(ret) <- dn
  ret
}

"is.rectangular" <- 
function(x)
{
  ( is(x, "character") || is(x, "numeric") || is(x, "complex") ||
   is(x, "logical") ||
   is(x, "factor") || is(x, "matrix") || is(x, "data.frame") ||
   (is(x, "array") && (length(dim(x)) == 2)) ||
   (is(x, "named") && (is.rectangular(x@.Data))) || is(x, "groupVec") ||
   is(x, "seriesVirtual") )
}

setMethod("rowIds", signature(x="ANY"),
          function(x)
          {
            dm <- dim(x)
            if(is.null(dm))
              return(names(x))
            dn <- dimnames(x)[[1]]
            if(length(dn) == dm[1])
              dn
          })

setReplaceMethod("rowIds", signature(x="ANY"),
                 function(x, value)
                 {
                   if((length(value) != 0) && (length(value) != numRows(x)))
                     stop("Incorrect length for rowIds")
                   dm <- dim(x)
                   if(is.null(dm))
                     names(x) <- value
                   else {
                     dn <- dimnames(x)
                     if(length(dn) < length(dm)) {
                       ## dimnames were NULL (or messed up)
                       if(is.null(value)) return(x)
                       dn <- vector("list", length(dm))
                     }
                     dn[1] <- list(value)
                     dimnames(x) <- dn
                   }
                   x
                 })

setMethod("colIds", signature(x="ANY"),
function(x)
{
  dm <- dim(x)
  if(!is.null(dm)) {
    dn <- dimnames(x)
    if(length(dn) < 2)
      return(NULL)
    else return(dn[[2]])
  }
  return(NULL)
})

setMethod("numCols", signature(x="ANY"),
function(x)
{
  nc <- ncol(x)
  if(is.null(nc) || is.na(nc))
    1
  else nc
})


setReplaceMethod("colIds", signature(x="ANY"),
                 function(x, value){
     if((length(value) != 0) && (length(value) != numCols(x)))
       stop("Incorrect length for colIds")
     
     dm <- dim(x)
     if(is.null(dm)) {
       x <- matrix(x)
     }
     if(is.null(value))
       value <- character(0)
     dn <- dimnames(x)
     if(length(dn) > 1)
       dn[[2]] <- value
     else dn <- list(character(0), value)
     dimnames(x) <- dn
     x
   })

setReplaceMethod("numCols", signature(x="data.frame"),
                 function(x, value){
                 ncol(x) <- value
                 return(x)
               })

setReplaceMethod("numCols", signature(x="ANY"),
  function(x, value){
    nc <- ncol(x)
    ids <- colIds(x)
    if(!is.null(nc))
      ncol(x) <- value
    else {
      if(value == 1L)
        return(x)
      else {
        x <- cbind(x, matrix(NA, length(x), value-1))
      }
    }
    if(!is.null(ids)) {
      numRows(ids) <- value
      colIds(x) <- ids
    }
    x
  })

setMethod("numRows", signature(x="ANY"),
   function(x)
   {
     nr <- nrow(x) 
     if(is.null(nr))
       length(x)
     else nr
   })

setReplaceMethod("numRows", signature(x="data.frame"),
                 function(x, value){
                 nrow(x) <- value
                 return(x)
               })

setReplaceMethod("numRows", signature(x="ANY"),
  function(x,  value){
    nr <- nrow(x)
    ids <- rowIds(x)
    if(!is.null(nr))
      nrow(x) <- value
    else {
      oldlen <- length(x)
      if(value == oldlen)
        return(x)
      else if (value <= oldlen)
        return(x[1:value])
      else{
        return(c(x, rep(NA, value-oldlen)))
      }
    }
    if(!is.null(ids)) {
      numRows(ids) <- value
      rowIds(x) <- ids
    }
    x
  })

setReplaceMethod("nrow",  signature("data.frame"),
function( x, value )
{
  olddim <- dim(x)
  if( value == olddim[1] )
    x
  else if( value <= olddim[1] )
    x[1:value,,drop=FALSE]
  else{
    x[value,] <- NA
    x
  }
})

setReplaceMethod("nrow",  signature("ANY"),
function( x, value )
## replacement function for nrow
## works for matrices, data.sheets and data.frames
{
  olddim <- dim(x)
  if( value == olddim[1] )
    x
  else if( value <= olddim[1] )
    x[1:value,,drop=FALSE]
  else
    rbind( x, matrix( NA, value - olddim[1], olddim[2] ))
})

setReplaceMethod("ncol", signature("data.frame"),
function(x, value)
  {
  olddim <- dim(x)
  if(value == olddim[2])
    x
  else if(value <= olddim[2])
    x[, 1:value, drop=FALSE]
  else{
    x[, (olddim[2] + 1):value] <- NA
    x
  }
  })

setReplaceMethod("ncol", signature("ANY"),
function(x, value)
{
  olddim <- dim(x)
  if(value == olddim[2])
    x
  else if(value <= olddim[2])
    x[, 1:value, drop = FALSE]
  else  cbind(x, matrix(NA, olddim[1], value - olddim[2]))
})
