##*****************************************************************************
## THIS IS NOT TO Roxygenize.

setClass("Grid",    	
         representation(
            dim = "integer",   
            dimNames = "character",
            levels = "list",
            index = "array"),
         validity = function(object) {
            if (length(object@dimNames) != object@dim) {
                return("incorrect length for 'dimNames'")
            }
            if (length(object@levels) != object@dim) {
                return("incorrect length for 'levels'")
            }
            if (!all(dim(object@index) == lapply(object@levels, length))) {
                return("'index' dimension conflicts with 'levels'")
            }
            if (any(sapply(object@levels, function(x) any(duplicated(x)) ))) {
                return("levels must be distinct'")
            }
       
            TRUE
         }
)

##*****************************************************************************
##' Create a Grid object.
##' 
##' @title Create a new \code{Grid} object
##'
##' @param levels A list with the levels of the variables. The length
##' of the list gives the spatial dimension. The \eqn{i}-th element
##' must be a numeric vector containing the distinct values of the
##' corresponding variable.
##'
##' @param nlevels Integer vector giving the number of levels by
##' dimension. This formal argument is used only when \code{levels} is
##' missing. In this case, the nodes for dimension \eqn{i} will be
##' regularly spaced between \code{0.0} and \code{1.0}.
##'
##' @param dimNames If the list \code{levels} does not have names,
##' then the character vector \code{dimNames} will be used.
##'
##' @param index An array with integer values and dimension
##' corresponding to the spatial dimensions. Each value represent the
##' number of the corresponding node in a "flat" representation
##' (data.frame or matrix), and will be used in tasks such as
##' interpolation.
##'
##' @return An object with S4 class \code{"Grid"}. If needed, this
##' object can be coerced to a data frame or to a matrix by using the
##' methods \code{as.data.frame} or \code{as.matrix}.
##'
##' @author Yves Deville
##'
##' @note If \code{index} is not provided, the vectors in the
##' \code{levels} list can be unsorted, and they will be
##' sorted. However, when a non-default \code{index} is provided, the
##' levels must be sorted. This rule is intended to reduce the risk of
##' error.
##'
##' @examples
##' myGD1 <- Grid(levels = list("X" = c(0, 1), "Y" = c(0.0, 0.5, 1.0)))
##' ## the same with less effort
##' myGD2 <- Grid(nlevels = c("X" = 2, "Y" = 3))
##' nlevels(myGD2)
##' levels(myGD2)
##' 
Grid <- function(levels,
                     nlevels = NULL,
                     dimNames = NULL,
                     index = NULL) {
  
  
  if (missing(levels)) {
    if ( (length(nlevels) == 1) && (!is.null(dimNames)) ) {
      nlevels <- rep(nlevels, length.out = length(dimnames))
      names(nlevels) <- dimNames
    }
    dimNames <- names(nlevels)
    levels <- list()
    for (i in 1L:length(nlevels)) {
      levels[[i]] <- seq(from = 0.0, to = 1.0, length.out = nlevels[i])
    }
    names(levels) <- names(nlevels)
  } else {
    if (!is.null(nlevels)){
      warning("'levels' is given so 'nlevels' is ignored")
    }
    if (!is.list(levels)) stop("'levels' must be a list")
    levels <- lapply(levels, as.numeric)
    if (!is.null(index) && any(sapply(levels, is.unsorted))) {
      stop("'index' can be provided only when levels are sorted")
    }
    levels <- lapply(levels, sort)
  }
  
  dim <- length(levels)
  nLevels <- sapply(levels, length)
  levChar <- lapply(levels, round_levels)
  levChar <- lapply(levChar, as.character)
      
  if (is.null(dimNames)) {
    if (!is.null(nm <- names(nLevels))) dimNames <- nm
    else {
      dimNames <- paste("X", 1L:dim, sep = "") 
      names(levels) <- dimNames
    }
  }
  names(levChar) <- dimNames
  
  n <- prod(nLevels)
  names(nLevels) <- NULL ## to avoid having names in the 'dim' attribute
  if (is.null(index)) index <- array(1L:n, dim = nLevels, dimnames = levChar)

  myGD <- new("Grid",
              dim = as.integer(dim),
              dimNames = as.character(dimNames),
              levels = levels,
              index = index)
  myGD
  
}

##*****************************************************************************
##' Random drawing of a \code{Grid} object
##'
##' The nodes are drawn by first choosing their number for each
##' spatial dimension using a Poisson distribution, and then by
##' drawing a sample of the uniform distribution with the suitable
##' size.
##' @title Random drawing of a \code{Grid} object
##'
##' @param dim The spatial dimension.
##'
##' @param nxMin Minimum number of nodes. Integer vector of length 1
##' or \code{dim}.example(Grid).
##'
##' @param nxE Expected number of nodes. Numeric vector of length 1
##' or \code{dim}.
##'
##' @param nlevels Integer vector giving the number of levels for each
##' spatial dimension.
##'
##' @param dimNames Names of the spatial dimensions.
##'
##' @return An object with S4 class \code{"Grid"}.
##'
##' @author Yves Deville
##'
##' @note This function is \emph{not} not for drawing random design
##' points from within the object.
##'
##' @seealso \code{\link{Grid-class}}.
##'
##' @examples
##' set.seed(1234)
##' ## specify dim: number of levels are random!
##' GD <- randGrid(dim = 4)
##' plot(GD)
##' ## specify number of levels (hence dim)
##' GD2 <- randGrid(nlevels = c(2, 3, 4))
##' plot(GD2)
##' 
randGrid <- function(dim = 1L + rpois(1L, lambda = 3),
                      nxMin = 1L, nxE = 4L,
                      nlevels, dimNames) {
   
   if (missing(nlevels)) {
     if (missing(dim))  {
       dim <- 1L + rpois(1L, lambda = 3)
     } 
     nlevels <- nxMin + rpois(dim, lambda = nxE - nxMin)
   } else {
     dim <- length(nlevels)
   }
   
   if (missing(dimNames)) {
      if (!is.null(nm <- names(nlevels))) dimNames <- nm
      else dimNames <- paste("X", 1L:dim, sep = "") 
   }
   
   levels <- list()
   for (i in 1L:dim) {
      levels[[i]] <- sort(runif(nlevels[i]))
   }
   names(levels) <- dimNames
   n <- prod(nlevels)
   
   levChar <- lapply(levels, round_levels)
   levChar <- lapply(levChar, as.character)
   
   index <- array(1L:n, dim = nlevels, dimnames = levChar)
   
   myGD <- new("Grid",
               dim = as.integer(dim),
               dimNames = as.character(dimNames),
               levels = levels,
               index = index)
   myGD
   
}

##*****************************************************************************
## Basic extraction methods
##*****************************************************************************

setMethod("dim", 
          signature = signature(x = "Grid"),
          definition = function(x){ x@dim })

setMethod("dimnames", 
          signature = signature(x = "Grid"),
          definition = function(x){ x@dimNames })

setMethod("dimnames<-", 
          signature = signature(x = "Grid"),
          definition = function(x, value){
              if (length(value) != x@dim) {
                  stop(sprintf("'value' has length %d but must have length %d",
                               length(value), x@dim))
              }
              x@dimNames[] <- names(x@levels) <- value
              L <- dimnames(x@index)
              names(L) <- value
              dimnames(x@index) <- L
              x
          }
          )

setMethod("levels", 
          signature = signature(x = "Grid"),
          definition = function(x){ x@levels })

setMethod("nlevels", 
          signature = signature(x = "Grid"),
          definition = function(x){
            sapply(x@levels, length)
          })

## XXX problem : define the 'index' function
## setMethod("index", 
##           signature = signature(object = "Grid"),
##           definition = function(object){ object@index })

##*****************************************************************************
## Coercion S4 methods
##*****************************************************************************
setMethod("as.data.frame", 
          signature = signature(x = "Grid"),
          definition = function(x, jitter = FALSE, ...){
             n <- prod(dim(x@index))
             lev <- x@levels
             names(lev) <- x@dimNames
             df <- expand.grid(lev)
             dni <- Matrix::invPerm(x@index)
             df <- df[dni, ]
             rownames(df) <- as.character(1L:n)
             if (jitter) {
                 dx <- lapply(lev, function(x) min(diff(x)) / 30)
                 for (j in 1:ncol(df)) {
                     if (is.finite(dx[[j]])) {
                         df[ , j] <-  df[ , j] + rnorm(n, mean = 0, sd = dx[[j]])
                     }
                 }
             }
             df
         }
)

setMethod("as.matrix", 
          signature = signature(x = "Grid"),
          definition = function(x, jitter = FALSE, ...){
              as.matrix(as.data.frame(x, jitter = jitter))
          })

##*****************************************************************************
## Classical methods
##*****************************************************************************          
setMethod("length", 
          signature = signature(x = "Grid"),
          definition = function(x){
              length(x@index)
          })


setMethod("plot", 
          signature = signature(x = "Grid"),
          definition = function(x, y, jitter = FALSE, ...){
              df <- as.data.frame(x, jitter = jitter)
              if (!missing(y)) {
                  df <- cbind(df, y)
              }
              pairs(df, ...)
          })

setMethod("show", 
          signature = signature(object = "Grid"),
          definition = function(object){
             cat("Grid Data object\n")
             cat(sprintf("  o dimension : %d\n", dim(object)))
             nx <- sapply(levels(object), length)
             nn <- sprintf("(%s)", dimnames(object))
             nn <- paste(nx, nn)
             nn <- paste(nn, collapse = ", ")
             cat(sprintf("  o dim names : %s\n",
                         paste(dimnames(object), collapse = ", ")))
             cat(sprintf("  o number of nodes : %s\n", nn))
             cat(sprintf("  o total number of nodes : %d\n", prod(nx)))
          })

##*****************************************************************************
## Permutation method
## The generalised transpose changes the order of the factors keeping the
## order of the observations. The transformation to a data frame is
## not efficient in term of computations, but is easily understood.
##*****************************************************************************          
## setMethod("aperm",
##           signature = signature(a = "Grid"),
##           definition = function(a, perm, ...){
##              indexNew <- aperm(a@index, perm = perm, ...)
##              dimNamesNew <- a@dimNames[perm]
##              levelsNew <- list()
##              for (j in 1L:dim(a)) {
##                  nm <- dimNamesNew[j]
##                  levelsNew[[nm]] <- a@levels[[perm[j]]] 
##              }
##              new("Grid",
##                  dim = dim(a),
##                  dimNames = dimNamesNew,
##                  levels = levelsNew,
##                  index = indexNew)      
##          })
setMethod("aperm",
          signature = signature(a = "Grid"),
          definition = function(a, perm, ...){
              if (!setequal(sort(as.integer(perm)), 1L:dim(a)))
                  stop("'perm' must be a valid permutation of ",
                       dim(a), " items")
              des <- as.data.frame(a)
              as.Grid(des[ , perm])
          })

##*****************************************************************************
## Coercion S3 methods
##*****************************************************************************
##' Coercion of objects into \code{Grid}.
##'
##' @title  Coercion of objects to \code{Grid}
##'
##' @param object An object to be coerced to \code{Grid}, typically a
##' data frame or a matrix.
##'
##' @param ... Not used yet. 
##'
##' @return An object with S4 class \code{"Grid"}.
as.Grid <- function(object, ...) {
   UseMethod("as.Grid")
}

##' Coercion to \code{Grid}.
##'
##' The dimensions of the \code{Grid} object are matched to the \eqn{d}
##' columns of the object given in \code{object}, in the same
##' order. Each column is coerced into a factor. Each combination
##' of the \eqn{d} different levels must be found exactly one in the
##' data frame or matrix, and the the row number of that combination is
##' stored in the \code{index} slot of the object.
##'
##' @title Coercion to \code{Grid}
##'
##' @aliases as.Grid.numeric as.Grid.matrix as.Grid.list
##'
##' @S3method as.Grid default
##' @S3method as.Grid data.frame
##'
##' @param object An object to be coerced into a \code{Grid} object.
##'
##' @param ... Not used yet.
##'
##' @return Object with S4 class \code{"Grid"}.
##'
##' @note A numeric object is coerced into a \code{Grid} object
##' with dimension \eqn{d =1}.
##'
##' @examples
##' set.seed(1234)
##' GDnum <- as.Grid(runif(8))
##' GDlist <- as.Grid(list("X" = runif(8), "Y" = runif(4)))
##' df <- expand.grid(X = runif(6), Y = runif(4))
##' GDdf <- as.Grid(df)
##' GDmat <- as.Grid(as.matrix(df))
##'
as.Grid.default <- function(object, ...) {
  as.Grid(as.data.frame(object), ...)
}

as.Grid.data.frame <- function(object, ...) {
   
   X <- lapply(object, as.factor)
   levels <- lapply(X, levels)
   levels <- lapply(levels, as.numeric)
   Table <- tapply(object[ , 1L], X, length)
   if (any(is.na(Table)) || any(Table != 1)) {
       stop("each combination of the levels must appear",
            " exactly once in the data frame 'object'")
   }
   n <- prod(sapply(levels, length))
   index <- tapply(1L:n, X, function(x) x)
   ## L <- dimnames(index)
   ## names(L) <- names(object)
   ## dimnames(index) <- L
   myGD <- new("Grid",
               dim = ncol(object),
               dimNames = names(object),
               levels = levels,
               index = index)
   myGD
   
}

as.Grid.list <- function(object, ...) {
  dim <- length(object)
  if (is.null(names(object))) {
    names(object) <- paste("X", 1L:dim, sep = "")
  }
  object <- lapply(object, as.numeric)
  object <- lapply(object, sort)
  nx <- sapply(object, length)
  index <- array(1L:prod(nx), dim = nx)
  myGD <- new("Grid",
               dim = dim,
               dimNames = names(object),
               levels = object,
               index = index)
   myGD
}
  
##*****************************************************************************
##' Apply a function to a \code{Grid} object.
##'
##' The function is applied to each vector combination of the levels by first
##' using the \code{as.matrix} coercion.
##'
##' @title Apply a function to a \code{Grid} object.
##' 
##' @param object An object with S4 class \code{"Grid"}
##'
##' @param fun The function to be applied.
##'
##' @param ... Further arguments to be passed to \code{fun}.
##'
##' @return A vector of values.
##'
##' @seealso \code{\link{array_Grid}} to reshape the result as an array.
##' 
##' @note The result values are given in the order specified by the
##' \code{Grid} object.
##' 
##' @examples
##' myGD <- Grid(levels = list("X" = c(0, 1),
##'                            "Y" = c(0.0, 0.5, 1.0),
##'                            "Z" = c(0.0, 0.2, 0.4, 1.0)))
##' ## gaussian function
##' apply_Grid(myGD, fun = function(x) exp(-x[1]^2 - 3 * x[2]^2 - 2 * x[3]^2))
##' 
##' 
apply_Grid <- function(object, fun, ...) {
    
  X <- as.matrix(object)
  fTest <- fun(X[1, ], ...)
  
  if (length(fTest) > 1L) {
    stop("'fun' returning a vector value is  not implemented yet")
  }

  Y <- apply(X = X, MARGIN = 1, FUN = fun, ...)
  
  ## is as.array is TRUE,
  if (FALSE) {
    dim(Y) <- nlevels(object)
    ## nms <- lapply(levels(object), as.character)
    dimnames(Y) <- dimnames(object@index)
  } 
  
  Y
  
}

##*****************************************************************************
##' Sample from/in a design object.
##'
##' From the object, \code{size} 'new' design points are sampled or
##' drawn. When \code{atSample} is \code{FALSE}, the range of each
##' column is computed or extracted from the object \code{x}, and then
##' independent drawings from the corresponding uniform distributions
##' are used.
##' 
##' @title Sample from/in a design object
##' 
##' @aliases sampleIn sampleIn,Grid-method sampleIn,data.frame-method
##' sampleIn,matrix-method
##' 
##'
##' @usage sampleIn(X, size = 1L, atSample = FALSE, ...)
##' 
##' @param X The object, typically a design matrix with class
##' \code{matrix} or an object with class \code{"Grid"}.
##'
##' @param size Number of sample points.
##'
##' @param atSample Logical. If \code{TRUE} the new sample points are
##' obtained using \code{sample}
##'
##' @param ... Other arguments to be passed to the \code{sample} function
##' when \code{atSample} is \code{TRUE}, typically \code{replace = TRUE}
##' is needed when the object \code{x} contains less than the target
##' number given in \code{size}.
##'
##' @return A matrix with \code{size} rows. When \code{atSample} is
##' \code{TRUE}, the matrix has an attribute named \code{index} which
##' gives the position of the sampled items in the original grid.
##'
##' @note When \code{x} is of class \code{"Grid"}, it may be the case
##' that the drawings with \code{atSample = FALSE} have no meaning since some
##' of the variables may be discrete variables and not continuous ones. 
##' 
##' @seealso The \code{\link{sample}} function.
##'
##' @examples
##' g <- Grid(levels = list("Temp" = c(400, 450, 500),
##'                         "Bore" = c(0, 600, 1200, 1800, 2400)))
##' X1 <- sampleIn(g, size = 4)
##' X2 <- sampleIn(g, size = 4, atSample = TRUE)
##' 
##' ## this must be zero
##' sum(abs(as.matrix(g)[attr(X2, "index"), ] - X2))
##'
##'
##'
##*****************************************************************************
## Sample
##*****************************************************************************
sampleIn <-  function(X, size = 1L, atSample = FALSE, ...) standardGeneric("sampleIn")
if (!isGeneric("sampleIn")) {
    setGeneric("sampleIn", sampleIn)
}

setMethod("sampleIn", 
          signature = signature(X = "matrix"),
          definition = function(X, size = 1L, atSample = FALSE, ...) {
              if (atSample) {
                  ind <- sample(x = 1L:nrow(X), size = size, ...)
                  XNew <- X[ind, , drop = FALSE]
                  attr(XNew, which = "index") <- ind
              } else {
                  XRange <- apply(X = X, MARGIN = 2, FUN = range)
                  XNew <- matrix(NA, nrow = size, ncol = ncol(X))
                  colnames(XNew) <- colnames(X)
                  for (j in 1:ncol(X)) {
                      XNew[ , j] <- runif(size, min = XRange[1, j],  max = XRange[2, j])
                  }
              }
              XNew   
          }
      )

setMethod("sampleIn", 
          signature = signature(X = "data.frame"),
          definition = function(X, size = 1L, atSample = FALSE, ...) {
              sampleIn(X = as.matrix(X), size = size, atSample = atSample, ...)
          }
          )

setMethod("sampleIn", 
          signature = signature(X = "Grid"),
          definition = function(X, size = 1L, atSample = FALSE, ...) {
              if (atSample) {
                  xMat <- as.matrix(X)
                  ind <- sample(x = 1L:nrow(xMat), size = size, ...)
                  XNew <- xMat[ind, , drop = FALSE]
                  attr(XNew, which = "index") <- ind
              } else {
                  XRange <- as.data.frame(lapply(levels(X), range))
                  XNew <- matrix(NA, nrow = size, ncol = dim(X))
                  colnames(XNew) <- dimnames(X)
                  for (j in 1L:dim(X)) {
                      XNew[ , j] <- runif(size, min = XRange[1, j],  max = XRange[2, j])
                  }
              }
              XNew   
          }
      )

##=======================================================================
##' Drop the dimensions with a unique level in a \code{Grid} object.
##'
##' @title Drop the dimensions with a unique level in a Grid
##' object
##' 
##' @param x An object with class \code{Grid} 
##'
##' @param ... Not used yet.
##'
##' @return An object of class \code{Grid} but with the fixed
##' dimensions removed.
##'
##' @seealso The usual \code{\link{drop}} function for arrays.
##' 
##' @examples
##' myGD <- Grid(nlevels = c(3, 1, 4, 1, 5))
##' myGD
##' plot(myGD)
##' myGD1 <- drop_Grid(myGD)
##  plot(myGD1)
drop_Grid <- function(x, ...) {
    nl <- nlevels(x) 
    fixed <- (nlevels(x) == 1L)
    if (any(fixed)) {
        levNew <- list()
        for (i in 1:dim(x)) {
            if (!fixed[i]) {
                nm <- x@dimNames[[i]]
                levNew[[nm]] <- x@levels[[i]]
            }
        }
        ## print(levNew)
        ind <- order(x@index)
        dim(ind) <- nl[!fixed]
        xNew <- new("Grid",
                    dim = x@dim - sum(fixed),   
                    dimNames = dimnames(x)[!fixed],
                    levels = levNew,
                    index =  ind)
        
        return(xNew)
        
    } else {
        return(x)
    }
    
}

##======================================================================
##' Find boundary points in a \code{Grid} object.
##'
##' When \code{type} is \code{"data.frame"}, the returned object is a
##' data frame containing only the boundary points. When \code{type}
##' is \code{"array"} the result is a logical array with its elements
##' in correspondence with the \code{index} slot of the
##' object. Finally, when \code{type} is \code{"index"}, the result is
##' an integer vector giving the indices of the boundary points in the
##' order of the nodes defined by the object. This order is the one
##' used by the \code{as.data.frame} coercion.
##' 
##' @title Find boundary points in a \code{Grid} object
##' 
##' @param x An object with class \code{"Grid"}.
##'
##' @param type The wanted type of result.
##'
##' @return An object describing the boundary points. See
##' \bold{Details}.
##'
##' @note Remind that when using the \code{plot} method for a
##' \code{Grid} object, some nodes are generally hidden because
##' several points have the same projection when they are shown in a
##' two-dimensional scatterplot.
##'
##' When one or more of the levels has length \eqn{\leq 2}{<= 2}, all
##' points of the grid are boundary points!
##'
##' @examples
##' ## define a Grid object
##' myGD <- Grid(nlevels = c(3, 4))
##' bd <- boundary_Grid(myGD, type = "index")
##'
##' ## use a different color for boundary points
##' cols <- rep("black", length(myGD))
##' cols[bd] <- "red"
##' plot(myGD, col = cols, pch = 16, cex = 2, main = "Boundary points")
##'
##' ## repeat this after a generalised transposition
##' myGD2 <- aperm(myGD, perm = c(2, 1))
##' bd2 <- boundary_Grid(myGD2, type = "index")
##' cols2 <- rep("black", length(myGD2))
##' cols2[bd2] <- "red"
##' plot(myGD2, col = cols2, pch = 16, cex = 2, main = "Boundary points")
##'
##' ## 3-dimensional
##' myGD3 <- Grid(nlevels = c("x" = 3, "y"= 4, "z" = 6))
##' bd3 <- boundary_Grid(myGD3, type = "index")
##' cols3 <- rep("black", length(myGD3))
##' cols3[bd3] <- "red"
##' plot(myGD3, jitter = TRUE, col = cols3, pch = 16, cex = 2, main = "Boundary points")
boundary_Grid <- function(x, type = c("data.frame", "array", "index")) {
    
    result <- match.arg(type)
    nl <- nlevels(x)
    n <- prod(nl)
    logA <- array(FALSE, dim = nlevels(x), dimnames = levels(x))
    
    for (i in 1:dim(x)) {
        si <- slice.index(x@index, i)
        logA[si == 1] <- TRUE
        logA[si == nl[i]] <- TRUE
    }
    
    if (result == "array") {
        return(logA)
    } else {
        dni <- as.vector(x@index)
        ## ind is the vector of number in the R-canonical order
        ind <- dni[logA]
        if (result == "index") return(ind)
        else if (result == "data.frame") {
            df <- as.data.frame(x)[ind, ]
            return(df)
        }
    }
    
}

##======================================================================
##' Subgrid by selection of levels in one dimension.
##'
##' @title Subgrid by selection of levels in one dimension.
##'
##' @param X An object with class \code{"Grid"}.
##'
##' @param subset An expression concerning \emph{one} of the dimension of \code{X}.
##' This can be compared to an expression in the \code{subset} method.
##' 
##' @param type The wanted type of result. When \code{type} is
##' \code{"index"}, the returned result is an integer vector
##' containing the indices of the subgrid nodes in the original
##' grid. When \code{type} is \code{"Grid"} the result is the subgrid
##' as an object of this class. Finally when \code{type} is
##' \code{"both"} a list is returned containing the two previous types
##' of result under the names \code{index} and \code{Grid}.
##' 
##' @param drop Logical. If \code{TRUE} and if only one level is selected
##' with \code{subset} then the corresponding dimension of the (flat)
##' \code{Grid} result will be dropped.
##'
##' @return A vector of integer indices, a \code{Grid} object or a list
##' embedding these two elements, depending on the value of \code{type}.
##'
##' @seealso The \code{\link{subset}} method.
##'
##' @note The new \code{Grid} returned (if any) uses consecutive
##' numbers for the nodes between \eqn{1} and the number of new nodes. So it
##' no longer keeps trace of the nodes position in \code{X}.  However
##' the surviving nodes are kept in the order that they had in
##' \code{X}.
##'
##' When a response vector or matrix say \code{Y} is related to the
##' original grid \code{X}, the index vector returned by
##' \code{subgrid} can be used to select the corresponding responses.
##' 
##' @examples
##' myGD <- Grid(levels = list("Temp" = c(250, 260, 280, 290, 300),
##'                            "Press" = c(3, 4, 6, 10)))
##' myGD2 <- subset_Grid(myGD, subset = Temp > 280) 
##' 
subset_Grid <- function(X, subset,
                    type = c("index", "Grid", "both"), drop = TRUE) {
  
    type <- match.arg(type)
    e <- substitute(subset)
    if (!is(X, "Grid")) {
        stop("'X' must be of class \"Grid\"")
    }
    av <- all.vars(e)
    if (length(av) > 1L) {
        stop("'subset' must concern only one variable")
    }
    j <- match(av, names(levels(X)))
    if (is.na(j)) {
        stop("variable '", av, "'does not match one of the dimension")
    }
    df <- as.data.frame(levels(X)[[j]])
    colnames(df) <- av
    ind <- eval(e, df)

    if (sum(ind) == 0) stop("empty subgrid")
    
    ind <- (1L:nlevels(X)[[j]])[ind]
    ## slice.index(X@index, j)
    ind2 <- X@index[slice.index(X@index, j) %in% ind]
    
    if (type == "index") return(ind2)
    else {
        Levels <- levels(X)
        Index <- Matrix::invPerm(order(ind2))
        
        if ((length(ind) == 1) && drop) {
            Levels[[j]] <- NULL
            dim(Index) <- sapply(Levels, length)
            X2 <- new("Grid",
                      dim = dim(X) - 1L,   
                      dimNames = dimnames(X)[-j],
                      levels = Levels,
                      index = Index)
                      
        } else {
            Levels[[j]] <- Levels[[j]][ind]
            dim(Index) <- sapply(Levels, length)
            X2 <- new("Grid",
                      dim = dim(X),
                      dimNames = dimnames(X),
                      levels = Levels,
                      index = Index)
        }
        if (type == "Grid") return(X2)
        else {
            return(list(index = ind2, Grid = X2))
        }
    }

}
