#' Iterative truncation procedure based on a bridge statistic.
#' @param x coda::mcmc sequence (will be cast to if necessary) to truncate transient
#' @param bridge bridge type to use: "brownian", student" or "loglik"
#' @param stat statistic to use for testing bridge: 
#' - if student bridge, "E","var","autocov","loglik_mean","loglik_extremum","extremum","ratio_extremum","ratio_loglik_extremum"
#' - if brownian bridge, "E","var","autocov","loglik_mean","loglik_extremum","extremum","ratio_extremum","ratio_loglik_extremum"
#' - if loglik bridge, "E","var","autocov","extremum","ratio_extremum"
#' @param param if "asymptotic" use asymptotic statistics, else if a list of 'N' and 'rho' use these parameters, if NULL estimate N and rho
#' @param trunc number of mcmc iterations to delete: if >=1, it is a constant number, if <1, a percentage of remaining batches
#' @param eps Target value for ratio of halfwidth to sample mean (for compatibility with heidel.diag)
#' @param pvalue significance level to use in iterative test
#' @seealso coda::heidel.diag
#' @references Heidelberger P and Welch PD. Simulation run length control in the presence of an initial transient. Opns Res., 31, 1109-44 (1983)
#' @examples
#' require(codadiags)
#' set.seed(123)
#' x = AR1()
#' print(bridgestat.diag(x))
#' y = add.transient(x)
#' print(bridgestat.diag(y,trunc=10))
bridgestat.diag <- function (x, bridge="student", stat="E", param="asymptotic", trunc=1, eps = 0.1, pvalue=0.3) {
    if (is.mcmc.list(x))
        return(lapply(x, bridgestat.diag, bridge, stat, trunc ,param, eps, pvalue))
    
    stat_str = paste(sep="",stat,".",bridge,"bridge")
    #if (is.character(bridge)) 
    bridge_str = paste(sep="",bridge,"bridge")
    bridge = eval(parse(text=bridge_str))
    #if (is.character(stat))
    stat = eval(parse(text=stat_str))
    
    
    eval.null.cdf <- function(x){
        if (is.null(param)) {
            rho = autocorr1(x)
            N = length(x)
            #parameters = paste("length =",N,", estimated auto-correlation =",rho)
            null.cdf = null.param.cdf(stat_str, N,rho)
        } else {
            if (is.character(param) && param == "asymptotic") {
                null.cdf = null.lim.cdf(stat_str)
                #parameters = "asymptotic distribution"
            } else {
                null.cdf = null.param.cdf(stat_str, param$N,param$rho)
                #parameters = paste("length =",param$N,", auto-correlation =",param$rho)
            }
        }
    }
    
    x <- as.mcmc(as.matrix(x))
    BS.mat0 <- matrix(0, ncol = 6, nrow = nvar(x))
    dimnames(BS.mat0) <- list(varnames(x),
                              c("stest", "start", "pvalue", "htest",
                                "mean", "halfwidth"))
    BS.mat <- BS.mat0
    for (j in 1:nvar(x)) {
        Y <- x[, j, drop = TRUE]    
        n1 <- length(Y)
        ## Schruben's-like test for convergence, applied sequentially on remaining batches
        # This is where things are really changing from heidel.diag
        while(niter(Y)>0) {
            n <- niter(Y)
            ybar <- mean(Y)
            b <- bridge(Y)
            STATISTIC <- stat(b)
            PVAL <- 1-eval.null.cdf(Y)(STATISTIC)
            if(PVAL > pvalue)
                break
            Y <- window(Y, start = start(Y)+ifelse(trunc<1,max(1,trunc*niter(Y)),trunc))
        }
        ## Recalculate S0 using section of chain that passed convergence test
        # for compatibility with coda
        S0ci <- spectrum0.ar(Y)$spec
        halfwidth <- 1.96 * sqrt(S0ci/n)
        passed <- !is.na(halfwidth) & (abs(halfwidth/ybar) <= eps)
        if (is.na(halfwidth)) {
            nstart <- NA
            passed <- NA
            halfwidth <- NA
            ybar <- NA
        }
        else {
            nstart <- start(Y)
        }
        converged <- TRUE # for compatibility with coda
        BS.mat[j, ] <- c(converged, nstart, PVAL, 
                         passed, ybar, halfwidth)
    }
    class(BS.mat) <- "bridgestat.diag"
    return(BS.mat)
}

# Print results of bridgstat.diag
# @param x bridgstat.diag call
# @param digits used to format values
# @param ... other arguments passed to print() calls
# @usage print(x, digits = 3, ...) 
# @aliases print,bridgestat.diag-method
"print.bridgestat.diag" <-
    function (x, digits = 3, ...) 
    {
        BS.title <- matrix(c("Stationarity", "test", "start", "iteration",
                             "p-value", "", 
                             "Halfwidth", "test", "Mean", "", "Halfwidth", ""),
                           nrow = 2)
        y <- matrix("", nrow = nrow(x), ncol = 6)
        for (j in 1:ncol(y)) {
            y[, j] <- format(x[, j], digits = digits)
        }
        y[, c(1, 4)] <- ifelse(x[, c(1, 4)], "passed", "failed")
        y <- rbind(BS.title, y)
        vnames <- if (is.null(rownames(x))) 
            paste("[,", 1:nrow(x), "]", sep = "")
        else rownames(x)
        dimnames(y) <- list(c("", "", vnames), rep("", 6))
        print.default(y[, 1:3], quote = FALSE, ...)
        print.default(y[, 4:6], quote = FALSE, ...)
        invisible(x)
    }





