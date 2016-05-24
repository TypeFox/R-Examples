##' @include guiContainer.R

##' Class for a gridded-layout container
setClass("gLayout",
         contains="guiContainer",
         prototype=prototype(new("guiContainer"))
         )

##' Constructor for grid layout container
##'
##' @export
glayout <- function(
                    homogeneous = FALSE, spacing = 10, container = NULL,      ... ,
                    toolkit=guiToolkit()){
  widget <- .glayout (toolkit,
                      homogeneous=homogeneous, spacing=spacing, container=container ,...
                      )
  obj <- new( 'gLayout',widget=widget,toolkit=toolkit) 
  return(obj)
}

##' generic for toolkit dispatch
##' @alias glayout
setGeneric( '.glayout' ,
           function(toolkit,
                    homogeneous = FALSE, spacing = 10, container = NULL,
                    ... )
           standardGeneric( '.glayout' ))

##' pass back item, list or matrix of items depending on dimension
setMethod("[",
          signature(x="gLayout"),
          function(x,i,j,...,drop=TRUE) {
            if(missing(drop)) drop <- TRUE
            if(missing(i))
              i <- seq_len(dim(x)[1])
            if(missing(j))
              j <- seq_len(dim(x)[2])

            if(length(i) == 1 && length(j) == 1)
              return(.leftBracket(x@widget, x@toolkit,i,j,...,drop=drop))

            ## a matrix or list
            out <- sapply(j, function(col) lapply(i, function(row) x[row, col]))
            if(is.matrix(out))
              return(out[,,drop=drop])
            else
              return(out)
             
          })

