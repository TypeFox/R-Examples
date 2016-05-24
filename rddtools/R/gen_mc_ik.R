#' Generate Monte Carlo simulations of Imbens and Kalyanaraman
#' 
#' Generate the simulations reported in Imbens and Kalyanaraman (2012)
#' @param n The size of sampel to generate
#' @param version The MC version of Imbens and Kalnayaraman (between 1 and 4).
#' @param sd The standard deviation of the error term.
#' @param output Whether to return a data-frame, or already a rdd_data
#' @param size The size of the effect, this depends on the specific version, defaults are as in ik: 0.04, NULL, 0.1, 0.1
#' @return An data frame with x and y variables. 
#' @export
#' @examples
#' mc1_dat <- gen_mc_ik()
#' MC1_rdd <- rdd_data(y=mc1_dat$y, x=mc1_dat$x, cutpoint=0)
#' 
#' ## Use np regression:
#' reg_nonpara <- rdd_reg_np(rdd_object=MC1_rdd)
#' reg_nonpara
#' 
#' # Represent the curves:
#' plotCu <- function(version=1, xlim=c(-0.1,0.1)){
#'   res <- gen_mc_ik(sd=0.0000001, n=1000, version=version)
#'   res <- res[order(res$x),]
#'   ylim <- range(subset(res, x>=min(xlim) & x<=max(xlim), 'y'))
#'   plot(res, type='l', xlim=xlim, ylim=ylim, main=paste('DGP', version))
#'   abline(v=0)
#'   xCut <- res[which(res$x==min(res$x[res$x>=0]))+c(0,-1),]
#'   points(xCut, col=2)
#' }
#' layout(matrix(1:4,2, byrow=TRUE))
#' plotCu(version=1)
#' plotCu(version=2)
#' plotCu(version=3)
#' plotCu(version=4)
#' layout(matrix(1))

gen_mc_ik <- function(n = 200, version = 1, sd = 0.1295, output = c("data.frame", "rdd_data"), size) {
    
    output <- match.arg(output)
    if (!version %in% c(1:4) | length(version) != 1) 
        stop("arg 'version' should be between 1 and 4")
    
    foo <- switch(version, `1` = gen_mc_ik_1, `2` = gen_mc_ik_2, `3` = gen_mc_ik_3, `4` = gen_mc_ik_4)
    if (missing(size)) {
        size <- switch(version, `1` = 0.04, `2` = 0, `3` = 0.1, `4` = 0.1)
    }
    res <- foo(n = n, sd = sd, size = size)
    if (output == "rdd_data") {
        res <- rdd_data(x = res$x, y = res$y, cutpoint = 0)
    }
    res
}


#################################### MC 1

gen_mc_ik_1 <- function(n = 200, sd = 0.1295, size = 0.04) {
    
    ## Regressor:
    Z <- rbeta(n, shape1 = 2, shape2 = 4, ncp = 0)
    X <- 2 * Z - 1
    error <- rnorm(n, sd = sd)
    
    ## Prepare variables:
    Y <- vector("numeric", length = n)
    ind_below <- X < 0
    X_low <- X[ind_below]
    X_up <- X[!ind_below]
    
    ## Compute Y variables:
    Y[ind_below] <- 0.48 + 1.27 * X_low + 7.18 * X_low^2 + 20.21 * X_low^3 + 21.54 * X_low^4 + 7.33 * X_low^5 + error[ind_below]
    Y[!ind_below] <- 0.48 + size + 0.84 * X_up - 3 * X_up^2 + 7.99 * X_up^3 - 9.01 * X_up^4 + 3.56 * X_up^5 + error[!ind_below]
    
    ## Result:
    res <- data.frame(x = X, y = Y)
    return(res)
}

#################################### MC 2

gen_mc_ik_2 <- function(n = 200, sd = 0.1295, size = 0) {
    
    # if(!missing(size) && !is.null(size)) warning('Argument 'size' ignored for gen_mc_ik_2') Regressor:
    Z <- rbeta(n, shape1 = 2, shape2 = 4, ncp = 0)
    X <- 2 * Z - 1
    error <- rnorm(n, sd = sd)
    
    ## Compute Y variables:
    Y <- ifelse(X < 0, 3 * X^2, 4 * X^2 + size) + error
    
    ## Result:
    res <- data.frame(x = X, y = Y)
    return(res)
}


#################################### MC 3

gen_mc_ik_3 <- function(n = 200, sd = 0.1295, size = 0.1) {
    
    ## Regressor:
    Z <- rbeta(n, shape1 = 2, shape2 = 4, ncp = 0)
    X <- 2 * Z - 1
    error <- rnorm(n, sd = sd)
    
    ## Compute Y variables:
    Y <- 0.42 + ifelse(X < 0, 0, size) + 0.84 * X - 3 * X^2 + 7.99 * X^3 - 9.01 * X^4 + 3.56 * X^5 + error
    
    ## Result:
    res <- data.frame(x = X, y = Y)
    return(res)
}

#################################### MC 4

gen_mc_ik_4 <- function(n = 200, sd = 0.1295, size = 0.1) {
    
    ## Regressor:
    Z <- rbeta(n, shape1 = 2, shape2 = 4, ncp = 0)
    X <- 2 * Z - 1
    error <- rnorm(n, sd = sd)
    
    ## Compute Y variables:
    Y <- 0.42 + ifelse(X < 0, 0, size) + 0.84 * X + 7.99 * X^3 - 9.01 * X^4 + 3.56 * X^5 + error
    
    ## Result:
    res <- data.frame(x = X, y = Y)
    return(res)
}


#################################### MC simple

gen_MC_simple <- function(n = 200, LATE = 0.3) {
    
    ## Regressor:
    x <- rnorm(n)
    D <- x >= 0
    y <- 0.8 + LATE * D + 0.3 * x + 0.1 * x * D + rnorm(n)
    rdd_data(x = x, y = y, cutpoint = 0)
    
} 
