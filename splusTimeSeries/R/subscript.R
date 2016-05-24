setMethod( "subscript2d" , signature(x="ANY"), 
          function(x, i, j)
{
  # Subscript function for rectangular objects
  lengthDim <- length(dim(x))
  if(lengthDim > 3)
    stop("subscript2d does not support arrays with 3+ dimensions")
  if(lengthDim == 2){
    if(!missing(i) && !missing(j))
      return(x[i, j, drop = FALSE])
    if(!missing(i))
      return(x[i,  , drop = FALSE])
    if(!missing(j))
      return(x[, j, drop = FALSE])
    return(x[,  , drop = FALSE])
  }

  # Rest is for atomic-like vectors
  if(!missing(j)) {
    if(mode(j) == "numeric") {
      j <- j[j != 0 & j != -1]
      if(!length(j))
        return(x[0])
      if(any(j != 1))
        stop("2nd subscript out of range")
      if(length(j) > 1)
        stop("2nd subscript too long")
    }
    else if((mode(j) == "logical") && (length(j) > 1))
      stop("2nd subscript too long")
    else stop("2nd vector subscript must be numeric or logical")
  }
  if(missing(i))
    return(x[])
  len <- length(x)
  if(mode(i) == "numeric") {
    if(!length(i))
      return(x[])
    if(any(i > len | i < -len ))
      stop("1st subscript out of range")
  }
  else if((mode(i) == "logical") && (length(i) > len))
    stop("1st subscript too long")
  else if(is(i, "character") && (is.null(names(x)) ||
                                 any(is.na(match(i, names(x))))))
    stop("non-matching 1st subscript")
  x[i]
})

setMethod( "subscript2d" , signature(x="data.frame"), 
          function(x, i, j)
{
  if(!missing(i) && !missing(j))
    return(x[i, j, drop = FALSE])
  if(!missing(i))
    return(x[i, , drop = FALSE])
  if(!missing(j))
    return(x[, j, drop = FALSE])
  x[,  , drop = FALSE]
})

setReplaceMethod("subscript2d", signature(x="ANY"),
function(x, i, j, value)
{
  ## Subscript replace function for rectangular objects
  ## default method works on atomic-like vectors
  if(!missing(j)) {
    if(mode(j) == "numeric") {
      j <- j[j != 0 & j != -1]
      if(!length(j))
        return(x)
      if(any(j != 1))
        stop("2nd subscript out of range")
      if(length(j) > 1)
        stop("2nd subscript too long")
    }
    else if((mode(j) == "logical") && (length(j) > 1))
      stop("2nd subscript too long")
    else stop("2nd vector subscript must be numeric or logical")
  }
  if(missing(i)) {
    x[] <- value
    return(x)
  }
  len <- length(x)
  if(mode(i) == "numeric") {
    i <- i[i != 0 & i >  - len]
    if(!length(i))
      return(x)
    if(any(i > len))
      stop("1st subscript out of range")
  }
  else if((mode(i) == "logical") && (length(i) > len))
    stop("1st subscript too long")
  else if(is(i, "character") && (is.null(names(x)) ||
                                 any(is.na(match(i, names(x))))))
    stop("non-matching 1st subscript")
  x[i] <- value
  x
})

setReplaceMethod("subscript2d", signature(x="data.frame"),
function(x, i, j, value)
{
  ## fix -- check args for out of bounds?
  ## fix -- allow value to be a sensible vector
  if(!missing(i) && !missing(j))
    x[i, j] <- value
  else if(missing(i) && !missing(j))
    x[, j] <- value
  else if(!missing(i) && missing(j))
    x[i,  ] <- value
  else x[] <- value
  x
})

setMethod("subscript2d", "matrix",
	  function(x,i,j)
	  {
	    if( !missing(i) && !missing(j))
	      return(x[i,j,drop=FALSE])
	    if( !missing(i))
	      return(x[i,,drop=FALSE])
	    if( !missing(j))
	      return(x[,j,drop=FALSE])
	    x[,,drop=FALSE]
	  })

setReplaceMethod("subscript2d", "matrix",
	  function(x,i,j,value)
	  {
	    if( !missing(i) && !missing(j))
	      x[i,j] <- value
	    else if( missing(i) && !missing(j))
	      x[,j] <- value
	    else if( !missing(i) && missing(j))
	      x[i,] <- value
	    else
	      x[] <- value
	    x
	  })

