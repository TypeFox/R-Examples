#' Testing for balanced covariates: equality of means with t-test
#' 
#' Tests equality of means by a t-test for each covariate, between the two full groups or around the discontinuity threshold 
#' 
#' @param object object of class rdd_data
#' @param bw a bandwidth
#' @param paired Argument of the \code{\link{t.test}} function: logical indicating whether you want paired t-tests.
#' @param var.equal Argument of the \code{\link{t.test}} function:  logical variable indicating whether to treat the two variances as being equal
#' @param p.adjust Whether to adjust the p-values for multiple testing. Uses the \code{\link{p.adjust}} function
#' @param \ldots currently not used
#' @return A data frame with, for each covariate, the mean on each size, the difference, t-stat and ts p-value. 
#' @author Matthieu Stigler <\email{Matthieu.Stigler@@gmail.com}>
#' @seealso \code{\link{covarTest_dis}} for the Kolmogorov-Smirnov test of equality of distribution
#' @examples
#' data(house)
#' 
#' ## Add randomly generated covariates
#' set.seed(123)
#' n_Lee <- nrow(house)
#' Z <- data.frame(z1 = rnorm(n_Lee, sd=2), 
#'                 z2 = rnorm(n_Lee, mean = ifelse(house<0, 5, 8)), 
#'                 z3 = sample(letters, size = n_Lee, replace = TRUE))
#' house_rdd_Z <- rdd_data(y = house$y, x = house$x, covar = Z, cutpoint = 0)
#' 
#' ## test for equality of means around cutoff:
#' covarTest_mean(house_rdd_Z, bw=0.3)
#' 
#' ## Can also use function covarTest_dis() for Kolmogorov-Smirnov test:
#' covarTest_dis(house_rdd_Z, bw=0.3)
#' 
#' ## covarTest_mean works also on regression outputs (bw will be taken from the model)
#' reg_nonpara <- rdd_reg_np(rdd_object=house_rdd_Z)
#' covarTest_mean(reg_nonpara)





#' @export
covarTest_mean <- function(object, bw = NULL, paired = FALSE, var.equal = FALSE, p.adjust = c("none", "holm", "BH", "BY", "hochberg", 
    "hommel", "bonferroni")) UseMethod("covarTest_mean")

#' @rdname covarTest_mean
#' @export
covarTest_mean.rdd_data <- function(object, bw = NULL, paired = FALSE, var.equal = FALSE, p.adjust = c("none", "holm", "BH", 
    "BY", "hochberg", "hommel", "bonferroni")) {
    
    cutpoint <- getCutpoint(object)
    covar <- getCovar(object)
    cutvar <- object$x
    
    covarTest_mean_low(covar = covar, cutvar = cutvar, cutpoint = cutpoint, bw = bw, paired = paired, var.equal = var.equal, 
        p.adjust = p.adjust)
    
}


#' @rdname covarTest_mean
#' @export
covarTest_mean.rdd_reg <- function(object, bw = NULL, paired = FALSE, var.equal = FALSE, p.adjust = c("none", "holm", "BH", "BY", 
    "hochberg", "hommel", "bonferroni")) {
    
    cutpoint <- getCutpoint(object)
    dat <- object$RDDslot$rdd_data
    covar <- getCovar(dat)
    cutvar <- dat$x
    if (is.null(bw)) 
        bw <- getBW(object)
    
    covarTest_mean_low(covar = covar, cutvar = cutvar, cutpoint = cutpoint, bw = bw, paired = paired, var.equal = var.equal, 
        p.adjust = p.adjust)
    
}


covarTest_mean_low <- function(covar, cutvar, cutpoint, bw = NULL, paired = FALSE, var.equal = FALSE, p.adjust = c("none", "holm", 
    "BH", "BY", "hochberg", "hommel", "bonferroni")) {
    
    p.adjust <- match.arg(p.adjust)
    
    ## subset
    if (!is.null(bw)) {
        isInH <- cutvar >= cutpoint - bw & cutvar <= cutpoint + bw
        covar <- covar[isInH, ]
        cutvar <- cutvar[isInH]
    }
    regime <- cutvar < cutpoint
    
    ## Split data
    covar_num <- sapply(covar, as.numeric)
    
    tests <- apply(covar_num, 2, function(x) t.test(x[regime], x[!regime], paired = paired, var.equal = var.equal))
    tests_vals <- sapply(tests, function(x) c(x[["estimate"]], diff(x[["estimate"]]), x[c("statistic", "p.value")]))
    
    ## Adjust p values if required:
    if (p.adjust != "none") 
        tests_vals["p.value", ] <- p.adjust(tests_vals["p.value", ], method = p.adjust)
    
    ## Print results
    res <- t(tests_vals)
    colnames(res)[3] <- "Difference"
    res
    
    
}




#' Testing for balanced covariates: equality of distribution
#' 
#' Tests equality of distribution with a Kolmogorov-Smirnov for each covariates, between the two full groups or around the discontinuity threshold 
#' 
#' @param object object of class rdd_data
#' @param bw a bandwidth
#' @param exact Argument of the \code{\link{ks.test}} function: NULL or a logical indicating whether an exact p-value should be computed.
#' @param p.adjust Whether to adjust the p-values for multiple testing. Uses the \code{\link{p.adjust}} function
#' @param \ldots currently not used
#' @return A data frame  with, for each covariate, the K-S statistic and its p-value. 
#' @author Matthieu Stigler <\email{Matthieu.Stigler@@gmail.com}>
#' @seealso \code{\link{covarTest_mean}} for the t-test of equality of means
#' @examples
#' data(house)
#' 
#' ## Add randomly generated covariates
#' set.seed(123)
#' n_Lee <- nrow(house)
#' Z <- data.frame(z1 = rnorm(n_Lee, sd=2), 
#'                 z2 = rnorm(n_Lee, mean = ifelse(house<0, 5, 8)), 
#'                 z3 = sample(letters, size = n_Lee, replace = TRUE))
#' house_rdd_Z <- rdd_data(y = house$y, x = house$x, covar = Z, cutpoint = 0)
#' 
#' ## Kolmogorov-Smirnov test of equality in distribution:
#' covarTest_dis(house_rdd_Z, bw=0.3)
#' 
#' ## Can also use function covarTest_dis() for a t-test for equality of means around cutoff:
#' covarTest_mean(house_rdd_Z, bw=0.3)
#' ## covarTest_dis works also on regression outputs (bw will be taken from the model)
#' reg_nonpara <- rdd_reg_np(rdd_object=house_rdd_Z)
#' covarTest_dis(reg_nonpara)

#' @export
covarTest_dis <- function(object, bw, exact = NULL, p.adjust = c("none", "holm", "BH", "BY", "hochberg", "hommel", "bonferroni")) UseMethod("covarTest_dis")

#' @rdname covarTest_dis
#' @export
covarTest_dis.rdd_data <- function(object, bw = NULL, exact = FALSE, p.adjust = c("none", "holm", "BH", "BY", "hochberg", "hommel", 
    "bonferroni")) {
    
    cutpoint <- getCutpoint(object)
    covar <- getCovar(object)
    cutvar <- object$x
    
    covarTest_dis_low(covar = covar, cutvar = cutvar, cutpoint = cutpoint, bw = bw, exact = exact, p.adjust = p.adjust)
    
}

#' @rdname covarTest_dis
#' @export
covarTest_dis.rdd_reg <- function(object, bw = NULL, exact = FALSE, p.adjust = c("none", "holm", "BH", "BY", "hochberg", "hommel", 
    "bonferroni")) {
    
    cutpoint <- getCutpoint(object)
    dat <- object$RDDslot$rdd_data
    covar <- getCovar(dat)
    cutvar <- dat$x
    if (is.null(bw)) 
        bw <- getBW(object)
    
    covarTest_dis_low(covar = covar, cutvar = cutvar, cutpoint = cutpoint, bw = bw, exact = exact, p.adjust = p.adjust)
    
}

covarTest_dis_low <- function(covar, cutvar, cutpoint, bw = NULL, exact = NULL, p.adjust = c("none", "holm", "BH", "BY", "hochberg", 
    "hommel", "bonferroni")) {
    
    p.adjust <- match.arg(p.adjust)
    
    ## subset
    if (!is.null(bw)) {
        isInH <- cutvar >= cutpoint - bw & cutvar <= cutpoint + bw
        covar <- covar[isInH, ]
        cutvar <- cutvar[isInH]
    }
    regime <- cutvar < cutpoint
    
    
    
    ## Split data
    covar_num <- sapply(covar, as.numeric)
    
    tests <- apply(covar_num, 2, function(x) ks.test(x[regime], x[!regime], exact = exact))
    tests_vals <- sapply(tests, function(x) x[c("statistic", "p.value")])
    
    ## Adjust p values if required:
    if (p.adjust != "none") 
        tests_vals["p.value", ] <- p.adjust(tests_vals["p.value", ], method = p.adjust)
    
    ## Print results
    res <- t(tests_vals)
    res
    
    
} 
