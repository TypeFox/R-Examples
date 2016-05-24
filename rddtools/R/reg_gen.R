#' General polynomial estimator of the regression discontinuity
#' 
#' Compute RDD estimate allowing a locally kernel weighted version of any estimation function
#' possibly on the range specified by bandwidth
#' @param rdd_object Object of class rdd_data created by \code{\link{rdd_data}}
#' @param covariates Formula to include covariates
#' @param order Order of the polynomial regression. 
#' @param bw A bandwidth to specify the subset on which the kernel weighted regression is estimated
#' @param weights Optional weights to pass to the lm function. Note this cannot be entered together with \code{bw}
#' @param slope Whether slopes should be different on left or right (separate), or the same.
#' @param covar.opt Options for the inclusion of covariates. Way to include covariates, either in the main regression (\code{include}) or as regressors of y in a first step (\code{residual}). 
#' @param fun The function to estimate the parameters
#' @param \ldots Further arguments passed to fun. See the example. 
#' @details This function allows the user to use a custom estimating function, instead of the traditional \code{lm()}. 
#' It is assumed that the custom funciton has following behaviour:
#' \enumerate{
#'   \item A formula interface, together with a \code{data} argument
#'   \item A \code{weight} argument
#'   \item A coef(summary(x)) returning a data-frame containing a column Estimate
#' }
#' Note that for the last requirement, this can be accomodated by writing a specific \code{\link{rdd_coef}} 
#' function for the class of the object returned by \code{fun}. 
#' @return An object of class rdd_reg_lm and class lm, with specific print and plot methods
#' @references TODO
#' @export rdd_gen_reg 
#' @examples
#' ## Step 0: prepare data
#' data(house)
#' house_rdd <- rdd_data(y=house$y, x=house$x, cutpoint=0)
#' 
#' ## Estimate a local probit:
#' house_rdd$y <- with(house_rdd, ifelse(y<quantile(y, 0.25), 0,1))
#' reg_bin_glm <- rdd_gen_reg(rdd_object=house_rdd, fun= glm, family=binomial(link='probit'))
#' print(reg_bin_glm)
#' summary(reg_bin_glm)
#'

rdd_gen_reg <- function(rdd_object, fun = glm, covariates = NULL, order = 1, bw = NULL, slope = c("separate", "same"), covar.opt = list(strategy = c("include", 
    "residual"), slope = c("same", "separate"), bw = NULL), weights, ...) {
    
    checkIsRDD(rdd_object)
    cutpoint <- getCutpoint(rdd_object)
    
    slope <- match.arg(slope)
    
    if (!missing(weights) & !is.null(bw)) 
        stop("Cannot give both 'bw' and 'weights'")
    
    
    ## Subsetting
    dat <- as.data.frame(rdd_object)
    
    if (!is.null(bw)) {
        weights <- ifelse(dat$x >= cutpoint - bw & dat$x <= cutpoint + bw, 1, 0)
    } else if (!missing(weights)) {
        weights <- weights
    } else {
        weights <- NULL
    }
    
    ## Construct data
    if (missing(weights)) 
        weights <- NULL
    dat_step1 <- model.matrix(rdd_object, covariates = covariates, order = order, bw = bw, slope = slope, covar.opt = covar.opt)
    
    ## Regression
    reg <- fun(y ~ ., data = dat_step1, weights = weights, ...)
    
    ## Return
    RDDslot <- list()
    RDDslot$rdd_data <- rdd_object
    reg$RDDslot <- RDDslot
    class(reg) <- c("rdd_reg_lm", "rdd_reg", class(reg))
    attr(reg, "PolyOrder") <- order
    attr(reg, "cutpoint") <- cutpoint
    attr(reg, "slope") <- slope
    attr(reg, "RDDcall") <- match.call()
    attr(reg, "bw") <- bw
    reg
}

rdd_gen_reg_old <- function(rdd_object, covariates = ".", bw = rdd_bw_ik(rdd_object), slope = c("separate", "same"), fun = glm, 
    ...) {
    
    slope <- match.arg(slope)
    checkIsRDD(rdd_object)
    if (!is.function(fun)) 
        stop("Arg 'fun' should be a function")
    cutpoint <- getCutpoint(rdd_object)
    
    ## Construct data
    dat <- as.data.frame(rdd_object)
    
    dat_step1 <- dat[, c("y", "x")]
    dat_step1$x <- dat_step1$x - cutpoint
    dat_step1$D <- ifelse(dat_step1$x >= 0, 1, 0)
    if (slope == "separate") {
        dat_step1$x_right <- dat_step1$x * dat_step1$D
    }
    
    ### Weights
    kernel_w <- Kernel_tri(dat_step1[, "x"], center = 0, bw = bw)
    
    ## Regression
    reg <- fun(y ~ ., data = dat_step1, weights = kernel_w, ...)
    
    ## Return
    class(reg) <- c("rdd_reg_gen", "rdd_reg", class(reg))
    attr(reg, "RDDcall") <- match.call()
    attr(reg, "cutpoint") <- cutpoint
    attr(reg, "bw") <- bw
    reg
} 
