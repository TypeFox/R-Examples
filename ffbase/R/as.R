#' Character Vectors
#'
#' The generic function \code{as.character} converts \code{ff} vectors to characters.\cr
#'
#' @rdname as.character.ff
#' @method as.character ff 
#' @example ../examples/as.R
#' @param x a \code{ff} vector
#' @param ... other parameters passed on to chunk
#' @return A factor \code{ff} vector of the same length of x.
#' @export 
#' @export as.character.ff
#' @seealso \code{\link[base]{as.character}}
as.character.ff <- function(x, ...){
	levs <- unique(x)[]
	levs <- levs[!is.na(levs)]
	res <- ff(vmode="integer", length = length(x), levels=as.character(levs))
	for (i in chunk(x, ...)){
    Log$chunk(i)
		res[i] <- as.character(x[i])		
	}
	res		
}

as.ff_matrix <- function(x, ...){
  UseMethod("as.ff_matrix")
}

as.ff_matrix.ffdf <- function(x, ...){
  result <- ff(NA, dim = dim(x), vmode = names(maxffmode(vmode(x)))[1])
  dimnames(result) <- dimnames(x)
  for(i in chunk(x)){
    Log$chunk(i)
    result[i, ] <- as.matrix(x[i, ])
  }
  result
}

# 
# as.big.matrix.ffdf <- function(x, separated=FALSE, backingfile=NULL, backingpath=NULL, descriptorfile=NULL, shared=TRUE){
#   maxvmode <- maxffmode(vmode(x))
#   type <- switch(maxvmode,
#                  boolean = "short", logical = "short", quad = "short", nibble = "short",
#                  byte = "short", ubyte = "short", short = "short", 
#                  ushort = "integer", integer = "integer",
#                  single = "double", double = "double",
#                  complex = "notyetimplemented",
#                  raw = "notyetimplemented",
#                  character = "notyetimplemented",
#                  "notyetimplemented"              
#   )
#   stopifnot(type %in% c("char","short","integer","double"))
#   containsfactors <- sapply(physical(x), is.factor.ff)
#   if(sum(containsfactors) > 0){
#     warning("factors codes are inserted in big.matrix, not the factor levels")
#     containsfactors <- names(containsfactors)[containsfactors == TRUE]
#     for(column in containsfactors){
#       levels(x[[column]]) <- NULL
#     }
#   }
#   ## Create the big matrix
#   y <- big.matrix(nrow=nrow(x), 
#                   ncol=ncol(x), type=type, 
#                   init=NULL, dimnames=dimnames(x), 
#                   separated=separated,
#                   backingfile=backingfile, 
#                   backingpath=backingpath,
#                   descriptorfile=descriptorfile, 
#                   shared=shared)
#   ## And fill it up
#   for(i in chunk(x)){
#     idx <- as.integer(as.hi(i))              
#     y[idx, ] <- as.matrix(ffbase:::ffdfget_columnwise(x, index=idx))
#   }
#   return(y)
# }


#' Trivial implementation, but very handy
#'
#' Coerce a ffdf object to an ffdf object.
#' @param x ffdf object
#' @param ... not used.
#' @export
#' @export as.ffdf.ffdf
#' @importFrom ff as.ffdf
as.ffdf.ffdf <- function(x, ...){
  x
}

#' As ram for an ffdf to get your ffdf as a data frame in RAM
#'
#' Load your ffdf object in RAM into a data.frame.
#' @param x an object of class ffdf
#' @param ... not used.
#' @return a data.frame in RAM
#' @export
#' @export as.ram.ffdf
#' @importFrom ff as.ram
as.ram.ffdf <- function(x, ...){
  x[, , drop=FALSE]
}
