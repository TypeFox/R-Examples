#' Scale a matrix-like object
#' @description Scales the columns of a numeric matrix.
#' @param X A numeric matrix.
#' @param type The scaling type. See details.
#' @details This function is similar to \code{scale} in \code{base} package, but it can also
#' normalize the columns of a matrix. Suppose \eqn{x} is one of the columns
#' in matrix \code{X}. The "\code{standardize}" is defined as
#' \deqn{standardize = (x - mean(x))/sd(x)}
#' and the "\code{normalize}" is defined as
#' \deqn{normalize = (x - min(x))/(max(x) - min(x)).}
#' The "\code{none}" is just to keep the original values. It is designed for the
#' \code{grpss} function in purpose.
#' @return A scaled numeric matrix.
#' @author  Debin Qiu, Jeongyoun Ahn
#' @importFrom stats sd
#' @examples x <- matrix(1:18, ncol = 3)
#' Scale(x)  # standardization
#' Scale(x, type = "normalize")  # normalization
#' Scale(x, type = "none")  # do nothing
#' @export
#' @keywords internal
Scale <- function(X,type = c("standardize","normalize","none")) {
  if (!is.numeric(X))
    stop("'X' must be a numeric matrix")
  type <- match.arg(type)
  if (type != "none") {
    uni.mean <- colMeans(X)
    uni.sd <- apply(X,2,sd)
    uni.min <- apply(X,2,min)
    uni.max <- apply(X,2,max)
    if (any(uni.max - uni.min == 0) || any(uni.sd == 0))
      stop("column(s) with constant values exist(s)")
  }
  as.matrix(switch(type, none = X,
                   standardize = sweep(sweep(X,2,uni.mean,"-"),2,uni.sd,"/"),
                   normalize = sweep(sweep(X,2,uni.min,"-"),2,uni.max - uni.min,"/"))
  )
}
