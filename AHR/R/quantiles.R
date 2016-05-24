dispStat <- function(x, tau, data, surv.fit.fun, surv.fit.param) {
    tmp <- surv.fit.fun(x, data, surv.fit.param)
    mapply(function(a, b) (a - (1 - tau))^2 / b, tmp$S, tmp$V / nrow(data), USE.NAMES=FALSE)
}

#' Estimate arbitrary quantiles of a survival distribution based on the (weighted) Kaplan-Meier
#'
#' @title wkmQuantile
#' @param tau number between 0 and 1 specifiying the quantile to estimate
#' @param formula an object of class '"formula"' specifying the conditional survival model
#' @param data data frame containing the variables in formula
#' @param conf.level confidence level (or NULL if no confidence interval should be calculated)
#' @param null.value true value of quantile or NULL if no p-value should be calculated
#' @param rr.subset logical vector defining subset of observations to use for response rate estimation (default: use all observations)
#' @return An object of class '"survQuantile"'
#' @references brookmeyer_confidence_1982
#' @export
#' @examples
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' Y <- pmin(T, C)
#' D <- T <= C
#' Z <- rep(c(0,1), c(100, 100))
#' wkmQuantile(0.5, Surv(Y, D) ~ strata(Z), data.frame(Y=Y, D=D, Z=Z))
wkmQuantile <- function(tau, formula, data, conf.level=0.95, null.value=NULL, rr.subset=rep(TRUE, nrow(data))) {
    if(!is.null(formula)) data <- parseFormula(formula, data, one.sample=TRUE)
    ## if is.null(formula) is TRUE assume that the variables in data are named V,Y,D,W
    
    fitQuantile(tau, data, conf.level, null.value, "wkm", list(var=TRUE, cov=FALSE, alpha=1, left.limit=FALSE, rr.subset=rr.subset))
}

fitQuantile <- function(tau, data, conf.level=0.95, null.value=NULL, method="wkm", surv.fit.param=list(cov=FALSE, alpha=1, left.limit=FALSE, rr.subset=rep(TRUE, nrow(data)))) {
    if(tau < 0 || tau > 1) stop("tau outside [0,1]")

    obj <- list()

    if(method == "wkm") obj$method <- paste("Quantile estimation based on (weighted) Kaplan-Meier estimator")
    else if(method == "aj") obj$method <- paste("Quantile estimation based on Aalen-Johansen estimator")
    else stop(paste("Unknown method:", method))

    ## get function object
    surv.fit.fun <- get(method)
    
    obj$alternative <- "two.sided"
    
    bisect <- function(S, tau, df) {
        n <- length(S)
        if(n == 1) S[1]
        else {
            m <- ceiling(n/2)
            if(df(S[m]) < tau) bisect(S[1:m], tau, df)
            else bisect(S[(m+1):n], tau, df)
        }
    }
    
    if(method == "aj") Ys <- sort(data$time)
    else Ys <- sort(data$Y)
    
    obj$estimate <- bisect(Ys, tau, function(x) surv.fit.fun(x, data, surv.fit.param)$S)
    if(tau == 0.5) names(obj$estimate) <- "median"
    else names(obj$estimate) <- paste0(100 * tau, "%-quantile")
    
    ## p-value
    if(!is.null(null.value)) {
        obj$null.value <- null.value
        if(tau == 0.5) names(obj$null.value) <- "median"
        else names(obj$null.value) <- paste0(100 * tau, "%-quantile")
        if(min(Ys) > null.value || null.value > max(Ys)) {
            warning("Cannot calculate p-value: null.value outside data range!")
            obj$p.value <- NA
        } else {
              Z <- dispStat(null.value, tau, data, surv.fit.fun, surv.fit.param)
              obj$statistic <- Z
              obj$p.value <- pchisq(Z, df=1, lower.tail=FALSE)
        }
    }
    
    ## confidence interval
    if(!is.null(conf.level)) {
        st.jumps <- dispStat(Ys, tau, data, surv.fit.fun, surv.fit.param)
        ca <- qchisq(conf.level, df=1)
        y <- which(st.jumps < ca)
        obj$conf.int <- c(Ys[min(y)], Ys[max(y)])
        attr(obj$conf.int, "conf.level") <- conf.level
    }

    class(obj) <- "survQuantile"
    obj
}


#' Compare quantiles of two independent samples (ratio or difference) based on (weighted-) Kaplan-Meier estimator
#'
#' @title wkmCompareQuantiles
#' @param tau number between 0 and 1 specifying the quantile
#' @param formula an object of class '"formula"' specifying the conditional survival model
#' @param data data frame containing the variables in formula
#' @param conf.level confidence level (or NULL if no confidence interval should be calculated)
#' @param null.value true value of quantile ratio or difference
#' @param method either '"ratio"' or '"difference"'
#' @param p.value if TRUE p.value will be calculated (requires null.value)
#' @return An object of class '"survQuantile"', i.e. a list containing the estimated quantiles, confidence interval and p.value (if p.value = TRUE)
#' @references su_nonparametric_1993
#' @export
#' @examples
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' Y <- pmin(T, C)
#' D <- T <= C
#' Z <- rep(c(0,1), c(100, 100)) # treatment indicator
#' wkmCompareQuantiles(0.5, Surv(Y, D) ~ Z, data.frame(Y=Y, D=D, Z=Z))
wkmCompareQuantiles <- function(tau, formula, data, conf.level=0.95, null.value=1, method="ratio", p.value=FALSE) {

    if(!is.null(formula)) data <- parseFormula(formula, data)
    ## if is.null(formula) is TRUE assume that the variables in data are named V,Y,D,W,Trt

    fitCompareQuantiles(tau, data, conf.level, null.value, method, p.value, "wkm", list(var=TRUE, cov=FALSE, alpha=1, left.limit=FALSE, rr.subset=rep(TRUE, nrow(data))))
}

fitCompareQuantiles <- function(tau, data, conf.level=0.95, null.value=1, method="ratio", p.value=FALSE, surv.fit.fun="wkm", surv.fit.param=list(rr.subset=rep(TRUE, nrow(data)), cov=FALSE, alpha=1, left.limit=FALSE)) {    
    grps <- levels(data$Trt)
    if(length(grps) != 2) stop("Need exactly two groups!")
    
    data1 <- data[data$Trt == grps[1], ]
    data2 <- data[data$Trt == grps[2], ]

    trt.sub <- data$Trt[surv.fit.param$rr.subset]

    surv.fit.param1 <- surv.fit.param
    surv.fit.param2 <- surv.fit.param

    surv.fit.param1$rr.subset <- surv.fit.param$rr.subset[trt.sub == grps[1]]
    surv.fit.param2$rr.subset <- surv.fit.param$rr.subset[trt.sub == grps[2]]
    
    tau1 <- fitQuantile(tau=tau, data=data1, conf.level=conf.level, null.value=NULL, method=surv.fit.fun, surv.fit.param=surv.fit.param1)
    tau2 <- fitQuantile(tau=tau, data=data2, conf.level=conf.level, null.value=NULL, method=surv.fit.fun, surv.fit.param=surv.fit.param2)

    ## get function object
    surv.fit.fun <- get(surv.fit.fun)
    
    n1 <- nrow(data1)
    n2 <- nrow(data2)
    
    obj <- list(estimate=c(tau1$estimate, tau2$estimate))
    
    ##obj$value <- switch(method, ratio = obj$quantile1 / obj$quantile2, difference = obj$quantile1 - obj$quantile2)
    
    if(tau == 0.5) {
        names(obj$estimate) <- c("median (group 1)", "median (group 2)")
        qstr <- "Median"
    } else {
          names(obj$estimate) <- c(paste0(100 * tau, "%-quantile (group 1)"),
                                   paste0(100 * tau, "%-quantile (group 2)"))
          qstr <- "Quantile"
      }

    if(!is.null(null.value)) obj$null.value <- null.value
    else {
        if(method == "ratio") obj$null.value <- 1
        else obj$null.value <- 0
    }
    
    if(method == "ratio") {
        obj$estimate[3] <- tau1$estimate / tau2$estimate
        names(obj$estimate)[3] <- paste(qstr, "ratio")
        names(obj$null.value) <- paste(qstr, "ratio")
    } else if(method == "difference") {
          obj$estimate[3] <- tau1$estimate - tau2$estimate
          names(obj$estimate)[3] <- paste(qstr, "difference")
          names(obj$null.value) <- paste(qstr, "difference")
      } else stop(paste("Unknown method:", method))
    
    ## events
    T1 <- data1$Y[data1$D == 1]
    T2 <- data2$Y[data2$D == 1]
    T1s <- sort(T1)
    T2s <- sort(T2)

    ## chi-square statistic
    stat <- function(S, V) (S - (1 - tau))^2 / V

    ## confidence interval
    if(!is.null(conf.level)) {
        x <- switch(method,
                    ratio = outer(T1s, T2s, "/"),
                    difference = outer(T1s, T2s, "-"))
        
        x1 <- surv.fit.fun(T1s, data1, surv.fit.param1)
        x2 <- surv.fit.fun(T2s, data2, surv.fit.param2)
        
        chi.sq1 <- stat(x1$S, x1$V / n1)
        chi.sq2 <- stat(x2$S, x2$V / n2)
        
        y <- x[outer(chi.sq1, chi.sq2, "+") <= qchisq(1 - conf.level, df=1, lower.tail=FALSE)]
        
        obj$conf.int <- c(min(y), max(y))
        attr(obj$conf.int, "conf.level") <- conf.level  
    }
    
    ## p-value
    if(p.value & !is.null(null.value)) {
        get.surv <- function(x, data, param) surv.fit.fun(times=x, data=data, param=param)
        G <- switch(method,
                    ratio = function(x) {
                        tmp1 <- get.surv(x, data1, surv.fit.param1)
                        tmp2 <- get.surv(null.value * x, data2, surv.fit.param2)
                        stat(tmp1$S, tmp1$V / n1) + stat(tmp2$S, tmp2$V / n2)
                    },
                    difference = function(x) {
                        tmp1 <- get.surv(x, data1, surv.fit.param1)
                        tmp2 <- get.surv(null.value + x, data2, surv.fit.param2)
                        stat(tmp1$S, tmp1$V / n1) + stat(tmp2$S, tmp2$V / n2)
                    })
        ## minimum should be close to tau1$quantile
        r <- optimize(G, c(tau1$conf.int[1], tau1$conf.int[2]))
        obj$p.value <- pchisq(r$objective, df=1, lower.tail=FALSE)
    }

    class(obj) <- "survQuantile"
    obj$method <- paste(qstr, method)
    obj$alternative <- "two.sided"
    
    obj
}

#' Print survQuantile object
#'
#' @title print.wkmQuantile
#' @param x an object of class '"survQuantile"'.
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to or from other methods.
#' @method print survQuantile
#' @export
print.survQuantile <- function(x, digits=3, ...) {
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep = "\n")
    cat("\n")

    out <- character()
    if(!is.null(x$statistic))
        out <- c(out, paste("z =",
                            format(round(x$statistic, 4))))
    if(!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = digits)
        out <- c(out, paste("p-value",
                            if(substr(fp, 1L, 1L) == "<") fp else paste("=",fp)))
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if(!is.null(x$alternative) && !is.null(x$p.value)) {
        cat("alternative hypothesis: ")
        if(!is.null(x$null.value)) {
            if(length(x$null.value) == 1L) {
                alt.char <-
                    switch(x$alternative,
                           two.sided = "not equal to",
                           less = "less than",
                           greater = "greater than")
                cat("true ", names(x$null.value), " is ", alt.char, " ",
                    x$null.value, "\n", sep = "")
            }
            else {
                cat(x$alternative, "\nnull values:\n", sep = "")
                print(x$null.value, ...)
            }
        }
        else cat(x$alternative, "\n", sep = "")
    }
    if(!is.null(x$estimate)) {
        cat("Estimates:\n")
        print(x$estimate, ...)
    }

    cat("\n")
    
    if(!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")),
            " percent confidence interval:\n", " ",
            paste(format(c(x$conf.int[1L], x$conf.int[2L])), collapse = " "),
            "\n\n", sep = "")
    }
    
    cat("\n")
    invisible(x)
}
