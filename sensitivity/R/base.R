# Base functions of the sensitivity package
# Gilles Pujol, 2006 - 2008
# Modified by Frank Weber (2016): support model functions
# returning a matrix or a 3-dimensional array.

ask <- function(x, ...)
  UseMethod("ask")


tell <- function(x, y = NULL, ...)
  UseMethod("tell")


response <- function(x, loop = FALSE, other_types_allowed = FALSE, ...) {
  id <- deparse(substitute(x))
  
  if (class(x$model) == "function") {
    if (loop) {
      n <- nrow(x$X)
      y <- sapply(1:n, function(i){
        x$model(x$X[i,], ...)
      }, simplify = "array")
      if(is.matrix(y)){
        if(typeof(y) == "list"){
          stop("The model function returns a list when applied to a row of x$X")
        } else{
          # This means the model function returned vectors, so y is a 
          # (m times n)-matrix (with unknown m). For better consistency with the 
          # default case (that y is a vector of length n), y is transposed to a 
          # (n times m)-matrix.
          y <- t(y)
        }
      } else if(is.array(y) && length(dim(y)) == 3){
        # This means the model function returned matrices of the same sizes.
        # Let m be the number of rows and z be the number of columns.
        # For better consistency with the other cases (y is a vector or a 
        # matrix), y is transformed to an array of dimensions c(n, m, z).
        
        # Change the order of the dimensions:
        y <- aperm(y, perm = c(3, 1, 2))
      } else if(!is.numeric(y)){
        stop("x$model returned an object that can't be handled")
      }
    } else {
      y <- x$model(x$X, ...)
    }
  } else if (TRUE %in% (paste("predict.", class(x$model), sep="") %in% methods(predict))) {
    y <- predict(x$model, x$X, ...)
  } else {
    stop("The model isn't a function or does not have a predict method")
  }
  
  if(other_types_allowed){
    if (!class(y)[1] %in% c("numeric", "matrix", "array") ||
        (is.array(y) && typeof(y) == "list")) {
      y <- as.numeric(y)
      warning("Conversion of the response to numeric")
    } else if(class(y) == "array" && length(dim(y)) > 3){
      stop("If the model returns an array, it must not have more ",
           "than 3 dimensions")
    }
  } else{
    if (class(y) != "numeric") {
      y <- as.numeric(y)
      warning("Conversion of the response to numeric")
    }
  }
  
  # Assign column names resp. dimnames if not existing:
  if(class(y) == "matrix" && 
     (is.null(colnames(y)) || any(nchar(colnames(y)) == 0))){
    colnames(y) <- paste0("ycol", 1:ncol(y))
  } else if(class(y) == "array"){
    if(is.null(dimnames(y))){
      dimnames(y) <- list(NULL, paste0("ycol", 1:dim(y)[2]), 
                          paste0("ydim3_", 1:dim(y)[3]))
    } else{
      if(is.null(dimnames(y)[[2]]) || any(nchar(dimnames(y)[[2]]) == 0)){
        dimnames(y)[[2]] <- paste0("ycol", 1:dim(y)[2])
      }
      if(is.null(dimnames(y)[[3]]) || any(nchar(dimnames(y)[[3]]) == 0)){
        dimnames(y)[[3]] <- paste0("ydim3_", 1:dim(y)[3])
      }
    } 
  }
  
  x$y <- y
  assign(id, x, parent.frame())
}

########################################
# function added in 2014

plotFG <- function(x) 
  UseMethod("plotFG")

