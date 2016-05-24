#' @encoding UTF-8
#' @title Unity-based normalization
#'
#' @description Normalizes as feature scaling \code{min - max}, or unity-based normalization typically used to bring the values into the range [0,1].
#'
#' @param x is a vector to be normalized.
#' @param method A string for the method used for normalization. Default is \code{method = "range"}, which brings the values into the range [0,1]. See details for other implemented methods.
#'
#' @return Normalized values in an object of the same class as \code{x}.
#' @details This approach may also be generalized to restrict the range of values to any arbitrary values \code{a}  and  \code{b}, using: \deqn{X' = a + \frac{(x - x_{min})(b - a)}{(x_{max} - x_{min})}}.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @seealso  \code{\link{svTransform}}, \code{\link{scale}}.
#'
#' @examples
#'
#' x <- sample(10)
#' normalize(x)
#'
#' # equivalently to
#' (x-min(x))/(max(x)-min(x))
#'
#' @keywords Rescaling
#' @keywords Transformation
#'
#' @export normalize
#' @docType methods
#' @rdname normalize-methods
#' @aliases normalize,numeric,character,ANY-method
`normalize`<-setClass("normalize", representation(x = "numeric", method="character"))

setGeneric("normalize", def=function(x, method = "range"){
  standardGeneric("normalize")
})

#' @rdname normalize-methods
setMethod(f="normalize", definition=function(x, method = "range"){
  method = .Match(arg = method, choices = c("range", "center", "Z-score", "z-score", "scale"))
  mat <- as.matrix(x)
  if(method=="range"){
  min_attr = apply(mat, 2, min)
  max_attr = apply(mat, 2, max)
  mat <- sweep(mat, 2, min_attr, FUN="-")
  ans = sweep(mat, 2,  max_attr-min_attr, "/")
  attr(ans, 'normalized:min') = min_attr
  attr(ans, 'normalized:max') = max_attr
  return(ans)
  }
  else if(method=="center"){
  }
  else if(method=="Z-score"||method=="z-score"||method=="scale"){
  }
  else if (!is.numeric(resul <- x))
    warning("Data not numeric, normalization not applicable")
  else stop("Unknown input method")
})
NULL


# check that we get mean of 0 and sd of 1
#colMeans(scaled.dat)  # faster version of apply(scaled.dat, 2, mean)
#apply(scaled.dat, 2, sd)


#`normalize.factor` <- function(x, range, domain=range(1:nlevels(x)), ...) {
#  width <- diff(range)
#  n <- length(levels(x)) - 1
#  range[1]  - 1/n + width * as.numeric(x) / n
#}

