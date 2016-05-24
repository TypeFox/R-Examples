#' RDD coefficient
#'
#' Function to access the RDD coefficient in the various regressions
#' @param object A RDD regression object
#' @param allInfo whether to return just the coefficients (allInfo=FALSE) or also the se/t stat/pval. 
#' @param allCo Whether to give only the RDD coefficient (allCo=FALSE) or all coefficients
#' @param \ldots Further arguments passed to/from specific methods
#' @return Either a numeric value of the RDD coefficient estimate, or a data frame with the estimate, 
#' its standard value, t test and p-value and 
#' @export


rdd_coef <- function(object, allInfo = FALSE, allCo = FALSE, ...) UseMethod("rdd_coef")

#' @rdname rdd_coef
#' @export
rdd_coef.default <- function(object, allInfo = FALSE, allCo = FALSE, ...) {
    res <- coef(summary(object))
    if (!allCo) 
        res <- res["D", , drop = FALSE]
    if (!allInfo) 
        res <- res[, "Estimate"]
    res
}

#' @rdname rdd_coef
#' @export
rdd_coef.rdd_reg_np <- function(object, allInfo = FALSE, allCo = FALSE, ...) {
    res <- object$coefMat
    if (!allCo) 
        res <- res["D", , drop = FALSE]
    if (!allInfo) 
        res <- res[, "Estimate"]
    res
}
 
