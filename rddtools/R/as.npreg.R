#' Convert an rdd_reg object to a  \code{npreg} object
#' 
#' Convert an rdd_object to a non-parametric regression \code{npreg} from package \code{np}
#' @param x Object of class \code{rdd_reg} created by \code{\link{rdd_reg_np}} or \code{\link{rdd_reg_lm}}
#' @param \ldots Further arguments passed to the \code{\link{npregbw}} or \code{\link{npreg}}
#' @details This function converts an rdd_reg object into an \code{npreg} object from package \code{np}
#' Note that the output won't be the same, since \code{npreg} does not offer a triangular kernel, but a Gaussian or Epanechinkov one. 
#' Another reason why estimates might differ slightly is that \code{npreg} implements a multivariate kernel, while rdd_reg 
#' proceeds as if the kernel was univariate. A simple solution to make the multivariate kernel similar to the  univariate one 
#' is to set the bandwidth for x and Dx to a large number, so that they converge towards a constant, and one obtains back the univariate kernel. 
#' @export
#' @return An object of class \code{npreg} or \code{npregbw}
#' @seealso \code{\link{as.lm}} which converts \code{rdd_reg} objects into \code{lm}.
#' @examples
#' # Estimate ususal rdd_reg:
#'  data(house)
#'  house_rdd <- rdd_data(y=house$y, x=house$x, cutpoint=0)
#'  reg_nonpara <- rdd_reg_np(rdd_object=house_rdd)
#' 
#' ## Convert to npreg:
#'  reg_nonpara_np <- as.npreg(reg_nonpara)
#'  reg_nonpara_np
#'  rdd_coef(reg_nonpara_np, allCo=TRUE, allInfo=TRUE)
#' 
#' ## Compare with result obtained with a Gaussian kernel:
#'  bw_lm <- dnorm(house_rdd$x, sd=rddtools:::getBW(reg_nonpara))
#'  reg_nonpara_gaus <- rdd_reg_lm(rdd_object=house_rdd, w=bw_lm)
#'  all.equal(rdd_coef(reg_nonpara_gaus),rdd_coef(reg_nonpara_np)) 


as.npregbw <- function(x, ...) {
    res <- as.npregbw_low(x = x, npreg = FALSE, ...)
    res
}

#' @rdname as.npregbw
#' @export
as.npreg <- function(x, ...) {
    res <- as.npregbw_low(x = x, npreg = TRUE, ...)
    res
}


as.npregbw_low <- function(x, npreg = FALSE, adjustik_bw = TRUE, ...) {
    
    dat <- getOriginalData(x)
    bw <- getBW(x)
    cutpoint <- getCutpoint(x)
    
    ## Specify inputs to npregbw:
    
    ## data:
    x <- dat$x
    dat_np <- data.frame(y = dat$y, x = x, D = ifelse(x >= cutpoint, 1, 0), Dx = ifelse(x >= cutpoint, x, 0))
    dataPoints <- data.frame(x = c(cutpoint, cutpoint), D = c(0, 1), Dx = c(0, cutpoint))
    
    ## bw:
    range.x <- range(dat$x, na.rm = TRUE, finite = TRUE)
    if (adjustik_bw) {
        ## & names(bw) =='h_opt'
        bw <- rdd_bw_ik(dat, kernel = "Normal")
    }
    bw_other <- 9999 * diff(range.x)
    bws <- c(bw, rep(bw_other, 2))
    
    
    ## start npregbw
    res <- np::npregbw(bws = bws, formula = y ~ x + D + Dx, data = dat_np, regtype = "ll", eval = dataPoints, bandwidth.compute = FALSE, 
        gradients = TRUE, ...)
    class(res) <- c("rdd_reg_npregbw", class(res))
    
    ## if npreg, return instead model_np <- npreg(bw_np, newdata=dataPoints, gradients=TRUE)
    if (npreg == TRUE) {
        
        # check if np is installed
        if (!requireNamespace("np", quietly = TRUE)) {
            stop("The package 'np' is needed for this function to work. Please install it.", call. = FALSE)
        }
        
        # require('np') requireNamespace('np', quietly = TRUE)
        options(np.messages = TRUE)  ## otherwise got warnings messages... probably because comes only if loaded!
        res <- np::npreg(res, newdata = dataPoints, gradients = TRUE, ...)
        class(res) <- c("rdd_reg_npreg", class(res))
    }
    
    attr(res, "RDDdf") <- dat_np
    attr(res, "cutpoint") <- cutpoint
    res
}


#' @export
rdd_coef.rdd_reg_npreg <- function(object, allInfo = FALSE, allCo = FALSE, ...) {
    
    co <- diff(object$mean)
    if (allInfo) {
        se <- sum(object$merr)
        zval <- co/se
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        res <- cbind(co, se, zval, pval)
        colnames(res) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
        rownames(res) <- "D"
    } else {
        res <- co
    }
    
    if (allCo) {
        cos <- c(object$mean[1], object$grad)
        ses <- c(object$merr[1], object$gerr)
        
        ## X_right:
        dataPoints_Xr <- data.frame(x = 0, D = 0, Dx = c(0, 1))
        Xr <- diff(predict(object, newdata = dataPoints_Xr))
        
        estimates <- c(cos[1], co, cos[2], Xr)
        
        if (allInfo) {
            zvals <- cos/ses
            pvals <- 2 * pnorm(abs(zvals), lower.tail = FALSE)
            res <- data.frame(Estimate = estimates, `Std. Error` = c(ses[1], se, ses[2:3]), `z value` = c(zvals[1], zval, zvals[2:3]), 
                `Pr(>|z|)` = c(pvals[1], pval, pvals[2:3]), check.names = FALSE)
            rownames(res) <- c("(Intercept)", "D", "x_left", "x_right")
        } else {
            res <- estimates
        }
    }
    
    res
} 
