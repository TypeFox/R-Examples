wavk.test <-
function(x, factor.length=c("user.defined", "adaptive.selection"), 
Window=round(0.1*length(x)), q=3/4, j=c(8:11), B=1000, 
H0=c("no trend", "linear"), method=c("boot", "asympt"), 
ar.order=NULL, BIC=TRUE, robust=TRUE, out=FALSE)
{
    ### Perform various checks.
    if (NCOL(x) > 1 | !is.numeric(x)) {
        stop("x is not a vector or univariate time series.")
    }
    n <- length(x)
    if (any(is.na(x))) {
        stop("x contains missing values.")
    }
    if (any(ls() == "factor.length")) {
        factor.length <- match.arg(factor.length)
    } else {
        factor.length <- "user.defined"
    }
    if (NCOL(q) > 1 | !is.numeric(q) | NROW(q) > 1) {
        stop("q is not a scalar.")
    }
    if (q >= 1 | q <= 0) {
        stop("q is out of range from 0 to 1.")
    }
    if (!is.vector(j) | !is.numeric(j)) {
        stop("j is not a numeric vector.")
    }
    if (factor.length == "user.defined") {
        kn <- Window[1]
    } else {
        kn <- length(x)*q^j
    }
    kn <- unique(sort(floor(kn)))
    kn <- kn[kn > 2 & kn < n]
    if (length(kn) == 0) {
        stop("set a proper window.")
    }
    if (factor.length == "adaptive.selection" & length(kn) < 3) {
        stop("number of possible windows is not enough for adaptive selection. Change parameters 'q' and/or 'j'.")
    }
    if (factor.length == "adaptive.selection") {
        method <- "boot"
    }
    B <- round(B)
    if (B <= 0) {
        stop("number of bootstrap samples B must be positive.")
    }
    H0 <- match.arg(H0)
    method <- match.arg(method)
    if (!is.null(ar.order) & (NCOL(ar.order) > 1 | !is.numeric(ar.order) | NROW(ar.order) > 1)) {
        stop("ar.order is not a scalar.")
    }
    if (!is.null(ar.order) && ar.order < 0) {
          stop("ar.order must be non-negative.")
    }

    ### Function.
    DNAME <- deparse(substitute(x))
    t <- c(1:n)/n
    result <- matrix(NA, length(kn), 2)
    res <- matrix(NA, 1, 2)
    if (is.null(ar.order)) {
        ar.order <- floor(10*log10(n))
    }

    if (H0 == "linear") {
        ALTERNATIVE <- "presence of a nonlinear trend"
        mod <- lm(x ~ t)
        linear <- mod$coefficients
        beta0 <- mod$coefficients[1]
        beta1 <- mod$coefficients[2]
        TS <- mod$residuals
    } else {
        H0 <- "no trend"
        ALTERNATIVE <- "presence of a trend" 
        TS <- x
    }
    TS <- as.vector(TS)

    if (ar.order > 0) {
        if (BIC) {
            bic <- rep(NA, ar.order)
            if (robust){
                for (i in 1:length(bic)) {
		        pheta_rob <- HVK(TS, ar.order=i)
		        tmp <- filter(TS, pheta_rob, sides=1)
		        et <- TS[(i+1):n] - tmp[i:(n-1)]
		        bic[i] <- (n-i)*log(var(et)) + i*log(n-i)
                }
                pheta <- HVK(TS, ar.order=which.min(bic))
            } else {
                for (i in 1:length(bic)) {
                    a <- ar(TS, aic=FALSE, order.max=i, demean=TRUE)
                    bic[i] <- n*log(a$var.pred) + i*log(n)
                }
                a <- ar(TS, aic=FALSE, order.max=which.min(bic), demean=TRUE)
                pheta <- a$ar
            }
        } else { #BIC==FALSE, use fixed ar.order>0
            if (robust){
                pheta <- HVK(TS, ar.order=ar.order)
            } else {
                a <- ar(TS, aic=FALSE, order.max=ar.order, demean=TRUE)
                pheta <- a$ar
            }
        }
        names(pheta) <- paste(rep("phi_", length(pheta)), c(1:length(pheta)), sep="")
    } else { #ar.order==0, no filtering
        pheta <- numeric(0)
    }
    
    if (length(pheta)>0) {
        tmp <- filter(x, pheta, sides=1)
        if (H0 == "linear") {
            tmp2 <- filter(mod$fitted.values, pheta, sides=1)
            Z <- (x[(length(pheta)+1):n] - tmp[length(pheta):(n-1)]) - (mod$fitted.values[(length(pheta)+1):n] - tmp2[length(pheta):(n-1)])
            ESTIMATE <- list(linear, pheta)
            names(ESTIMATE) <- c("linear_trend_coefficients", "AR_coefficients")
        } else {
            Z <- (x[(length(pheta)+1):n] - tmp[length(pheta):(n-1)])
            ESTIMATE <- list(pheta)
            names(ESTIMATE) <- c("AR_coefficients")
        }
    } else {
        Z <- TS
        if (H0 == "linear") {
            ESTIMATE <- list(linear, pheta)
            names(ESTIMATE) <- c("linear_trend_coefficients", "AR_coefficients")
        } else {
            ESTIMATE <- list(pheta)
            names(ESTIMATE) <- c("AR_coefficients")
        }
    }

    Z <- Z - mean(Z)
    sigma <- sqrt(sum(diff(Z)^2)/(2*(length(Z)-1)))
    if (method == "asympt") {
        for (i in 1:length(kn)) {
            tmp <- WAVK(Z, kn[i])
            result[i,] <- c(tmp$Tns, tmp$p.value)
        }
        STATISTIC <- result[1,1]
        P.VALUE <- result[1,2]
        PARAMETER <- kn[1]
        names(PARAMETER) <- "user-defined window"
    } else {
        boot <- array(data=rnorm(n*B), c(n,B)) * sigma
        s <- array(data=NA, c(length(kn), B))
        for (i in 1:length(kn)) {
            s[i,] <- apply(boot, 2, function(x) WAVK(x, kn[i])$Tns)
            result[i,1] <- WAVK(Z, kn[i])$Tns
            crit <- sum(result[i,1] < s[i,])/B
            if (crit < 0.5) {
                result[i,2] <- 2*crit
            } else {
                result[i,2] <- 2*(1-crit)
            }
        }
        if (length(kn) < 3) {
            STATISTIC <- result[1,1]
            P.VALUE <- result[1,2]
            PARAMETER <- kn[1]
            names(PARAMETER) <- "user-defined window"
        } else {
            distance <- rep(NA, length(kn)-1)
            for (i in 1:length(distance)) {
                distance[i] <- dist(rbind(sort(s[i,]), sort(s[i+1,])))
            }
            kn_opt <- kn[which.min(distance)]
            res[1,] <- result[which.min(distance),]
            STATISTIC <- res[1,1]
            P.VALUE <- res[1,2]
            PARAMETER <- kn_opt
            names(PARAMETER) <- "adaptively selected window"
            if (out) {
                tmp <- names(ESTIMATE)
                ESTIMATE <- c(ESTIMATE, NA)
                names(ESTIMATE) <- c(tmp, "all_considered_windows")
                ESTIMATE[[length(ESTIMATE)]] <- cbind(kn, result)
                dimnames(ESTIMATE[[length(ESTIMATE)]]) <- list(rep("", length(kn)), c("Window", "WAVK-statistic", "p-value"))
            }
        }
    }
    METHOD <- "Trend test by Wang, Akritas and Van Keilegom"
    names(STATISTIC) <- "WAVK test statistic"
    if (out) {
        structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, p.value = P.VALUE, alternative = ALTERNATIVE, parameter = PARAMETER, estimate = ESTIMATE), class = "htest") 
    } else {
        structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, p.value = P.VALUE, alternative = ALTERNATIVE, parameter = PARAMETER), class = "htest") 
    }
}
