#' Parametric polynomial estimator of the regression discontinuity
#' 
#' Compute a parametric polynomial regression of the ATE, 
#' possibly on the range specified by bandwidth
#' @param rdd_object Object of class rdd_data created by \code{\link{rdd_data}}
#' @param covariates TODO
#' @param bw A bandwidth to specify the subset on which the parametric regression is estimated
#' @param inference Type of inference to conduct: non-parametric one (\code{np}) or standard (\code{lm}). See details. 
#' @param slope Whether slopes should be different on left or right (separate), or the same.
#' @param covar.opt Options for the inclusion of covariates. Way to include covariates, either in the main regression (\code{include}) or as regressors of y in a first step (\code{residual}). 
#' @return An object of class rdd_reg_np and class lm, with specific print and plot methods
#' @seealso \code{\link{rdd_bw_ik}} Bandwidth selection using the plug-in bandwidth of Imbens and Kalyanaraman (2012)
#' @references TODO
#' @export rdd_reg_np
#' @examples
#' ## Step 0: prepare data
#' data(house)
#' house_rdd <- rdd_data(y=house$y, x=house$x, cutpoint=0)
#' ## Step 2: regression
#' # Simple polynomial of order 1:
#' reg_nonpara <- rdd_reg_np(rdd_object=house_rdd)
#' print(reg_nonpara)
#' plot(reg_nonpara)


rdd_reg_np <- function(rdd_object, covariates = NULL, bw = rdd_bw_ik(rdd_object), slope = c("separate", "same"), inference = c("np", 
    "lm"), covar.opt = list(slope = c("same", "separate"), bw = NULL)) {
    
    slope <- match.arg(slope)
    inference <- match.arg(inference)
    checkIsRDD(rdd_object)
    cutpoint <- getCutpoint(rdd_object)
    
    if (!is.null(covariates)) 
        warning("covariates not fully implemented for non-para reg")
    
    ## Construct data
    if ("strategy" %in% names(covar.opt)) 
        warning("Arg 'strategy' should not be used for ")
    covar.opt$strategy <- "include"
    dat <- as.data.frame(rdd_object)
    dat_step1 <- model.matrix(rdd_object, covariates = covariates, order = 1, bw = bw, slope = slope, covar.opt = covar.opt)
    
    
    ### Weights
    kernel_w <- Kernel_tri(dat_step1[, "x"], center = 0, bw = bw)
    
    ## Regression
    reg <- lm(y ~ ., data = dat_step1, weights = kernel_w)
    coefD <- coef(reg)["D"]
    
    ## Non-para inference:
    if (inference == "np") {
        var <- var_estim(x = dat$x, y = dat$y, point = cutpoint, bw = bw, eachSide = TRUE)
        dens <- dens_estim(x = dat$x, point = cutpoint, bw = bw, eachSide = TRUE)
        
        const <- 4.8/(nrow(dat) * bw)
        all <- const * sum(var)/dens
        se <- sqrt(all)
        tval <- coefD/se
        pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
        coefmat <- matrix(c(coefD, se, tval, pval), nrow = 1, dimnames = list("D", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
    } else {
        coefmat <- coef(summary(reg))  #['D', , drop=FALSE]
    }
    
    ## Return
    res <- list()
    RDDslot <- list()
    RDDslot$rdd_data <- rdd_object
    RDDslot$model <- reg
    res$coefficients <- coef(reg)["D"]
    res$coefMat <- coefmat
    res$residuals <- residuals(reg)
    res$fitted <- fitted(reg)
    res$RDDslot <- RDDslot
    
    class(res) <- c("rdd_reg_np", "rdd_reg", "lm")
    attr(res, "RDDcall") <- match.call()
    attr(res, "cutpoint") <- cutpoint
    attr(res, "bw") <- bw
    res
}


#' @export 
print.rdd_reg_np <- function(x, signif.stars = getOption("show.signif.stars"), ...) {
    
    RDDcall <- attr(x, "RDDcall")
    bw <- getBW(x)
    cutpoint <- getCutpoint(x)
    x_var <- getOriginalX(x)
    
    n_left <- sum(x_var >= cutpoint - bw & x_var < cutpoint)
    n_right <- sum(x_var >= cutpoint & x_var <= cutpoint + bw)
    
    cat("### RDD regression: nonparametric local linear###\n")
    cat("\tBandwidth: ", bw, "\n")
    cat("\tNumber of obs: ", sum(n_left + n_right), " (left: ", n_left, ", right: ", n_right, ")\n", sep = "")
    
    cat("\n\tCoefficient:\n")
    
    printCoefmat(rdd_coef(x, allInfo = TRUE), signif.stars = signif.stars)
    
}


#' @export
summary.rdd_reg_np <- function(object, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), 
    ...) {
    
    x <- object
    bw <- getBW(x)
    cutpoint <- getCutpoint(x)
    x_var <- getOriginalX(x)
    
    ## compute numbers left/right:
    n_left <- sum(x_var >= cutpoint - bw & x_var < cutpoint)
    n_right <- sum(x_var >= cutpoint & x_var <= cutpoint + bw)
    
    ## compute residual summary:
    res_quant <- quantile(residuals(x))
    names(res_quant) <- c("Min", "1Q", "Median", "3Q", "Max")
    
    ## compute R^2
    r.squared <- summary(x$RDDslot$model)$r.squared
    
    ## Extend the rdd_reg_no output with new computaations:
    
    object$r.squared <- r.squared
    object$res_quant <- res_quant
    object$n_obs <- list(n_left = n_left, n_right = n_right, total = n_left + n_right)
    
    class(object) <- c("summary.rdd_reg_np", class(object))
    object
}


#' @export
print.summary.rdd_reg_np <- function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), 
    ...) {
    
    bw <- getBW(x)
    
    cat("### RDD regression: nonparametric local linear###\n")
    cat("\tBandwidth: ", bw, "\n")
    cat("\tNumber of obs: ", x$n_obs$total, " (left: ", x$n_obs$n_left, ", right: ", x$n_obs$n_right, ")\n", sep = "")
    
    cat("\n\tWeighted Residuals:\n")
    print(zapsmall(x$res_quant, digits + 1))
    
    
    cat("\n\tCoefficient:\n")
    
    printCoefmat(rdd_coef(x, allInfo = TRUE), signif.stars = signif.stars)
    
    cat("\n\tLocal R squared:", formatC(x$r.squared, digits = digits), "\n")
    
}


#' @export
plot.rdd_reg_np <- function(x, binwidth, chart = c("locpoly", "np"), ...) {
    
    chart <- match.arg(chart)
    cutpoint <- getCutpoint(x)
    bw <- getBW(x)
    if (missing(binwidth)) 
        binwidth <- bw/5  # binwidth!=bandwidth
    
    ## data
    dat <- getOriginalData(x, classRDD = FALSE)
    
    ## Use locpoly:
    dat_left <- subset(dat, x < cutpoint)
    dat_right <- subset(dat, x >= cutpoint)
    
    if (chart == "locpoly") {
        llp_left <- locpoly(x = dat_left$x, y = dat_left$y, bandwidth = bw)
        llp_right <- locpoly(x = dat_right$x, y = dat_right$y, bandwidth = bw)
        
        ## Use np:
    } else {
        np_reg_left <- np::npreg(np::npregbw(y ~ x, data = dat_left, regtype = "ll", ckertype = "epanechnikov", bandwidth.compute = FALSE, 
            bws = bw))
        
        np_reg_right <- np::npreg(np::npregbw(y ~ x, data = dat_right, regtype = "ll", ckertype = "epanechnikov", bandwidth.compute = FALSE, 
            bws = bw))
        newDat_left <- data.frame(x = seq(min(dat_left$x), cutpoint - 0.001, by = 0.01))
        newDat_right <- data.frame(x = seq(cutpoint, max(dat_right$x), by = 0.01))
        pred_left <- predict(np_reg_left, newdata = newDat_left, se.fit = TRUE)
        pred_right <- predict(np_reg_right, newdata = newDat_right, se.fit = TRUE)
    }
    ## plot
    plotBin(dat$x, dat$y, h = binwidth, ...)
    if (chart == "locpoly") {
        lines(llp_left$x, llp_left$y)
        lines(llp_right$x, llp_right$y)
    } else {
        lines(newDat_left$x, pred_left$fit, col = 1)
        lines(newDat_left$x, pred_left$fit + 2 * pred_left$se.fit, col = 2, lty = 2)
        lines(newDat_left$x, pred_left$fit - 2 * pred_left$se.fit, col = 2, lty = 2)
        
        lines(newDat_right$x, pred_right$fit, col = 1)
        lines(newDat_right$x, pred_right$fit + 2 * pred_right$se.fit, col = 2, lty = 2)
        lines(newDat_right$x, pred_right$fit - 2 * pred_right$se.fit, col = 2, lty = 2)
    }
}

#' @export 
vcov.rdd_reg_np <- function(object, ...) {
    
    infType <- infType(object)
    if (infType == "np") {
        warning("No vcov() available when rdd_reg_np() was called with infType='np'")
        res <- NULL
    } else {
        res <- vcov(object$RDDslot$model)
    }
    res
} 
