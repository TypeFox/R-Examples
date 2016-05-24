## rFormula
## Methods: formula, model.frame, model.matrix

#' Model formula for Rchoice models
#' 
#' Two kind of variables are used in models with individual heterogenetiy: the typical
#' variables that enter in the latent process and those variables that enter in the random
#' parameter (Hierarchical Model). \code{rFormula} deal with this type of models using
#' suitable methods to extract the elements of the model.
#' 
#' @param object a formula form the \code{rFormula} function, for the \code{model.matrix} method, a \code{rFormula} object,
#' @param formula a \code{rFormula} object,
#' @param data a \code{data.frame},
#' @param lhs see \code{\link[Formula]{Formula}},
#' @param rhs see \code{\link[Formula]{Formula}},
#' @param ... further arguments.
rFormula <- function(object){
  UseMethod("rFormula")
}

#' @rdname rFormula
#' @export
is.rFormula <- function(object){
  inherits(object, "rFormula")
}

rFormula.formula <- function(object){
  if (!inherits(object, "Formula")) object <- Formula(object)
  class(object) <- c("rFormula", "Formula", "formula")
  object
}

rFormula <- function(object){
  stopifnot(inherits(object, "formula"))
  if (!inherits(object, "Formula")) object <- Formula(object)
  if (!inherits(object, "rFormula")) class(object) <- c("rFormula", class(object))
  object
}

as.Formula.rFormula <- function(x, ...){
  class(x) <- class(x)[-1]
  x
}

#' @rdname rFormula
#' @export
model.frame.rFormula <- function(formula, data, ..., lhs = NULL, rhs = NULL){
  if (is.null(rhs)) rhs <- 1:(length(formula)[2])
  if (is.null(lhs)) lhs <- ifelse(length(formula)[1] > 0, 1, 0)
  index <- attr(data, "index")
  mf <- model.frame(as.Formula(formula), as.data.frame(data), ..., rhs = rhs)
  if(!is.null(index)) rownames(index) <- rownames(mf)
  index <- index[rownames(mf), ]
  index <- data.frame(lapply(index , function(x) x[drop = TRUE]), row.names = rownames(index))
  structure(mf,
            index = index,
            class = c("pdata.frame", class(mf)))
}

## has.intercept
has.intercept <- function(object, ...) {
  UseMethod("has.intercept")
}

has.intercept.default <- function(object, ...) {
  has.intercept(formula(object), ...)
}

has.intercept.formula <- function(object, ...) {
  attr(terms(object), "intercept") == 1L
}

has.intercept.Formula <- function(object, rhs = NULL, ...) {
  if (is.null(rhs)) rhs <- 1:length(attr(object, "rhs"))
  sapply(rhs, function(x) has.intercept(formula(object, lhs = 0, rhs = x)))
}

has.intercept.rFormula <- function(object, ...){
  attr(object, "class") <- "Formula"
  has.int <- has.intercept(object,...)
  has.int
}


## model matrix

#'@rdname rFormula
#'@export
model.matrix.rFormula <- function(object, data, rhs = NULL, ...){
  index <- attr(data, "index")
  if (is.null(rhs)) rhs <- 1
  
  if (rhs == 1){
    formula <- formula(object, rhs = 1, lhs = 0) # Normal covariates
    X <- model.matrix(formula, data)
  }
  if (rhs == 2){
    for.ind.esp <- formula(object, rhs = 2, lhs = 0)
    has.int <- has.intercept(for.ind.esp, rhs = 2)
    if (has.int) for.ind.esp <- update(for.ind.esp, ~ . -1)
    if(length(index) != 0L){
      id <- index[[1]]
      indata <- data[!duplicated(id), ]
      X <- model.matrix(for.ind.esp, indata)
    } else X <- model.matrix(for.ind.esp, data)
  }
  X
}  
 
is.hierarchical <- function(object) {
  ifelse(length(object)[2] == 2, TRUE, FALSE)
}
