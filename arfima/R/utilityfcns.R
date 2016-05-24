removeMode <- function(object, num) {
    if (!object$weeded) 
        stop("Please weed the object first before manually removing a mode.")
    m <- length(object$modes)
    if (num > m) 
        stop("incorrect mode specification")
    modes <- object$modes
    modes[[num]] <- NULL
    object$modes <- modes
    object <- weed(object, eps2 = 0)
    object
}

BIC <- function(object, ...) {
    UseMethod("BIC")
}

print.predarfima <- function(x, digits = max(6, getOption("digits") - 3), ...) {
    
    n <- length(x$z)
    
    m <- length(x) - 7
    
    seed <- x$seed
    limiting <- x$limiting
    bootpred <- x$bootpred
    B <- x$B
    predint <- x$predint
    
    ret <- vector("list", m)
    nam <- paste("Mode", 1:m)
    for (i in 1:m) {
        
        
        ans <- rbind(as.numeric(x[[i]]$Forecast), x[[i]]$exactSD, x[[i]]$limitSD)
        n.ahead <- dim(ans)[2]
        rownames(ans) <- c("Forecasts", "Exact SD", if (limiting) "Limiting SD" else NULL)
        colnames(ans) <- 1:n.ahead
        if (bootpred) {
            ans2 <- rbind(x[[i]]$uppernp, x[[i]]$meanvalnp, x[[i]]$lowernp)
            intt <- paste(round(100 * predint), "%", sep = "")
            namer1 <- c(paste("Upper", intt), "Prediction (Mean)", paste("Lower", intt))
            rownames(ans2) <- namer1
            colnames(ans2) <- 1:n.ahead
        } else ans2 <- NULL
        ret[[i]] <- list(`Forecasts and SDs` = ans, `Bootstrap Replicates` = B, `Bootstrap Predictions and Intervals` = ans2)
        if (B == 0) 
            ret[[i]][[2]] <- ret[[i]][[3]] <- NULL
    }
    names(ret) <- nam
    
    print(ret, digits = digits, ...)
    
    invisible(x)
}

print.arfima <- function(x, digits = max(6, getOption("digits") - 3), ...) {
    dmean <- x$dmean
    itmean <- x$itmean
    
    if (!x$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- x
    cat("Number of modes:", length(ans$modes), "\n")
    coeffs <- NULL
    nm <- namer <- NULL
    r <- ans$r
    s <- ans$s
    b <- ans$b
    
    
    period <- ans$period
    for (i in 1:length(ans$modes)) {
        phi <- ans$modes[[i]]$phi
        theta <- ans$modes[[i]]$theta
        phiseas <- ans$modes[[i]]$phiseas
        thetaseas <- ans$modes[[i]]$thetaseas
        phip <- ans$modes[[i]]$phip
        thetap <- ans$modes[[i]]$thetap
        phiseasp <- ans$modes[[i]]$phiseasp
        thetaseasp <- ans$modes[[i]]$thetaseasp
        
        dfrac <- ans$modes[[i]]$dfrac
        dfs <- ans$modes[[i]]$dfs
        H <- ans$modes[[i]]$H
        Hs <- ans$modes[[i]]$Hs
        alpha <- ans$modes[[i]]$alpha
        alphas <- ans$modes[[i]]$alphas
        
        delta <- ans$modes[[i]]$delta
        omega <- ans$modes[[i]]$omega
        logl <- ans$modes[[i]]$loglik
        sigma2 <- ans$modes[[i]]$sigma2
        muHat <- ans$modes[[i]]$muHat
        
        se <- ans$modes[[i]]$se
        if (i == 1) {
            nm <- getnames(ans)
            nm <- c(nm, paste("logl", sep = ""), paste("sigma^2", sep = ""))
            counter <- 0
            if (length(phi) > 1) {
                nm <- c(nm, paste("phi_p(", 1:length(phi), ")", sep = ""))
                counter <- counter + length(phi)
            }
            if (length(theta) > 1) {
                nm <- c(nm, paste("theta_p(", 1:length(theta), ")", sep = ""))
                counter <- counter + length(theta)
            }
            if (length(phiseas) > 1) {
                nm <- c(nm, paste("phi_p.", period, "(", 1:length(phiseas), ")", sep = ""))
                counter <- counter + length(phiseas)
            }
            if (length(thetaseas) > 1) {
                nm <- c(nm, paste("theta_p.", period, "(", 1:length(thetaseas), ")", sep = ""))
                counter <- counter + length(thetaseas)
            }
        }
        
        
        coeff <- c(phi, theta, phiseas, thetaseas, dfrac, H, alpha, dfs, Hs, alphas, muHat,
                   omega, delta, logl, sigma2)
        
        if (length(phi) > 1) 
            coeff <- c(coeff, phip)
        if (length(theta) > 1) 
            coeff <- c(coeff, thetap)
        if (length(phiseas) > 1) 
            coeff <- c(coeff, phiseasp)
        if (length(thetaseas) > 1) 
            coeff <- c(coeff, thetaseasp)
        if (x$getHess) {
            ses <- c(se, -Inf, -Inf, if (counter > 0) rep(-Inf, counter) else numeric(0))
            
            coeff <- rbind(coeff, ses)
        } else coeff <- matrix(coeff, nrow = 1)
        
        coeff <- signif(coeff, digits = digits)
        
        coeffs <- rbind(coeffs, coeff)
        star <- if (i %in% ans$index) 
            "*" else ""
        nam <- if (x$getHess) 
            c(paste("Coef.", i, star, ":", sep = ""), paste("SE.", i, star, ":", sep = "")) else paste("Coef.", i, star, ":", sep = "")
        namer <- c(namer, nam)
    }
    cat("\nCall:", deparse(x$call, width.cutoff = 75L), "", sep = "\n")
    cat("Coefficients for fits:\n", sep = "")
    
    colnames(coeffs) <- nm
    rownames(coeffs) <- namer
    coeffs <- t(coeffs)
    indss <- which(!is.finite(coeffs) & !is.na(coeffs), arr.ind = TRUE)
    inds2 <- which(is.na(coeffs), arr.ind = TRUE)
    coeffs2 <- coeffs
    pos <- which(is.finite(coeffs) & !is.na(coeffs) & (coeffs >= 0), arr.ind = TRUE)
    coeffs2[indss] <- ""
    coeffs2[inds2] <- "NA"
    coeffs2[pos] <- paste(" ", coeffs[pos], sep = "")
    print.default(coeffs2, print.gap = 2, quote = FALSE)
    cat("Starred fits are close to invertibility/stationarity boundaries\n")
    if (ans$numcut > 0) 
        cat("Wall detection on: ", ans$numcut, " mode(s) were believed to be spurious and were cut\n", 
            sep = "")
    if (any(is.na(coeffs))) 
        cat("NAs come from singular Hessians\n")
    invisible(x)
}

fitted.arfima <- function(object, ...) {
    
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    
    
    m <- length(ans$modes)
    
    ret <- vector("list", m)
    
    for (i in 1:m) {
        ret[[i]] <- ans$modes[[i]]$fitted
    }
    
    names(ret) <- paste("Mode", 1:m, sep = "")
    ret
}

residuals.arfima <- resid.arfima <- function(object, reg = FALSE, ...) {
    
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    
    
    xreg <- object$xreg
    if (reg && is.null(xreg)) 
        stop("Regression residuals requested, but no regression performed.")
    
    
    m <- length(ans$modes)
    
    ret <- vector("list", m)
    
    for (i in 1:m) {
        cmode <- ans$modes[[i]]
        if (!reg) 
            ret[[i]] <- cmode$residuals else ret[[i]] <- cmode$regResiduals
    }
    
    names(ret) <- paste("Mode", 1:m, sep = "")
    ret
}


summary.arfima <- function(object, digits = max(4, getOption("digits") - 3), ...) {
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    cov <- TRUE
    cor <- TRUE
    FInfo <- TRUE
    cf <- coef(ans)
    vcovs <- vcov(ans, type = "b", summ = TRUE)
    m <- length(ans$modes)
    dmean <- ans$dmean
    cff <- vector("list", m)
    if (FInfo) 
        finfo <- vector("list", m) else finfo <- NULL
    if (cov && cor) 
        corrs <- vcov(ans, cor = cor, type = "b", summ = TRUE) else if (cov) 
        corrs <- vcovs else corrs <- NULL
    logl <- logLik(ans)
    aics <- AIC(logl)
    bics <- AIC(logl, k = length(object$z))
    ans1 <- list()
    sigmas <- rep(0, m)
    tseflag <- TRUE
    for (i in 1:m) {
        ses <- sqrt(diag(vcovs[[i]]$observed))
        if ((!is.logical(dmean)) || (!dmean)) 
            ses <- c(NA, ses)
        if (!vcovs$warnH && !vcovs$warnX) 
            tses <- c(sqrt(diag(vcovs[[i]]$expected)), NA) else {
            tses <- NULL
            tseflag <- FALSE
        }
        if (!is.null(finfo)) 
            finfo[[i]] <- solve(vcovs[[i]]$expected * ans$n)
        zval <- cf[i, ]/ses
        ps <- 2 * pnorm(-abs(zval))
        cff[[i]] <- cbind(Estimate = cf[i, ], `Std. Error` = ses, `Th. Std. Err.` = tses, 
            `z-value` = zval, `Pr(>|z|)` = ps)
        sigmas[i] <- ans$modes[[i]]$sigma2
    }
    ans1$tseflag <- tseflag
    if (cor) 
        ans1$corrs <- corrs else ans1$covs <- corrs
    ans1$sigma2 <- sigmas
    ans1$logl <- logl
    ans1$coef <- cff
    ans1$aics <- aics
    ans1$bics <- bics
    ans1$finfo <- finfo
    ans1$m <- m
    ans1$call <- ans$call
    class(ans1) <- "summary.arfima"
    ans1
}

print.summary.arfima <- function(x, digits = max(6, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), 
    ...) {
    cat("\nCall:\n ", deparse(x$call, width.cutoff = 75L), "", sep = "\n")
    
    for (i in 1:x$m) {
        cat("\nMode", i, "Coefficients:\n")
        printCoefmat(x$coef[[i]], digits = digits, signif.stars = signif.stars, ...)
        cat(paste("sigma^2 estimated as ", format(x$sigma2[i], digits = digits), "; Log-likelihood = ", 
            format(x$logl[i], digits = digits), "; AIC = ", format(x$aics[i], digits = digits), 
            "; BIC = ", format(x$bics[i], digits = digits), "\n", sep = ""))
        if (!is.null(x$corrs)) {
            cat("\nNumerical Correlations of Coefficients:\n")
            correl <- format(round(x$corrs[[i]]$observed, 2), nsmall = 2, digits = digits)
            
            print(correl, quote = FALSE)
            cat("\nTheoretical Correlations of Coefficients:\n")
            correl <- format(round(x$corrs[[i]]$expected, 2), nsmall = 2, digits = digits)
            print(correl, quote = FALSE)
        } else if (!is.null(x$covs)) {
            cat("\nNumerical Covariances of Coefficients:\n")
            correl <- format(round(x$covs[[i]]$observed, 2), nsmall = 2, digits = digits)
            
            print(correl, quote = FALSE)
            cat("\nTheoretical Covariances of Coefficients:\n")
            correl <- format(round(x$covs[[i]]$expected, 2), nsmall = 2, digits = digits)
            print(correl, quote = FALSE)
        }
        if (!is.null(x$finfo)) {
            ff <- format(round(x$finfo[[i]], 2), nsmall = 2, digits = digits)
            cat("\nExpected Fisher Information Matrix of Coefficients:\n")
            print(ff, quote = FALSE)
        }
    }
    if (!x$tseflag) 
        cat("\nTheoretical standard errors not reported since there is either FGN or (dynamic) regression in the model\n")
    invisible(x)
}

InfoMatrix <- function(object, digits = max(6, getOption("digits") - 3), tapprox = FALSE, 
    ...) {
    if (class(object) != "arfima") 
        stop("only arfima objects can be passed to this function")
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    
    m <- length(ans$modes)
    
    ret <- vector("list", m)
    ifdf <- length(ans$modes[[1]]$dfrac) > 0
    ifdfs <- length(ans$modes[[1]]$dfs) > 0
    period <- ans$period
    for (i in 1:m) {
        cmode <- ans$modes[[i]]
        if (length(cmode$H) > 0 || length(cmode$Hs) > 0 || length(cmode$alpha) > 0 || length(cmode$alphas) > 
            0) 
            warning("Please note FGN and HD do not have known information matrix components: \n\n\t\treporting matrix without those components")
        ret[[i]] <- iARFIMA(phi = cmode$phi, theta = cmode$theta, phiseas = cmode$phiseas, 
            thetaseas = cmode$thetaseas, dfrac = ifdf, dfs = ifdfs, period = period, exact = !tapprox)
    }
    names(ret) <- paste("Mode", 1:m)
    ret
}



logLik.arfima <- function(object, ...) {
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    logl <- rep(0, length(ans$modes))
    for (i in 1:length(ans$modes)) logl[i] <- ans$modes[[i]]$loglik
    
    len <- ncol(coef(ans))
    
    if (is.numeric(ans$dmean) && !ans$itmean) 
        len <- len - 1
    
    # sigma^2 is the +1 in the below.
    attr(logl, "df") <- len + 1
    attr(logl, "nobs") <- attr(logl, "nall") <- ans$n
    class(logl) <- "logLik"
    logl
}

# confint.arfima

coef.arfima <- function(object, tpacf = FALSE, digits = max(4, getOption("digits") - 3), 
    ...) {
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    
    dmean <- ans$dmean
    itmean <- ans$itmean
    m <- length(ans$modes)
    s <- object$s
    r <- object$r
    b <- object$b
    coeffs <- NULL
    period <- ans$period
    for (i in 1:m) {
        phi <- ans$modes[[i]]$phi
        theta <- ans$modes[[i]]$theta
        phiseas <- ans$modes[[i]]$phiseas
        thetaseas <- ans$modes[[i]]$thetaseas
        phip <- ans$modes[[i]]$phip
        thetap <- ans$modes[[i]]$thetap
        phiseasp <- ans$modes[[i]]$phiseasp
        thetaseasp <- ans$modes[[i]]$thetaseasp
        
        dfrac <- ans$modes[[i]]$dfrac
        dfs <- ans$modes[[i]]$dfs
        H <- ans$modes[[i]]$H
        Hs <- ans$modes[[i]]$Hs
        alpha <- ans$modes[[i]]$alpha
        alphas <- ans$modes[[i]]$alphas
        
        delta <- ans$modes[[i]]$delta
        omega <- ans$modes[[i]]$omega
        
        muHat <- ans$modes[[i]]$muHat
        
        xreglist <- ans$modes[[i]]$xreglist
        if (i == 1) {
            nm <- NULL
            
            if (tpacf) {
                if (length(phip) > 0) 
                  nm <- c(nm, paste("phi_p(", 1:length(phi), ")", sep = ""))
                
                if (length(thetap) > 0) 
                  nm <- c(nm, paste("theta_p(", 1:length(theta), ")", sep = ""))
                
                if (length(phiseasp) > 0) 
                  nm <- c(nm, paste("phi_p.", period, "(", 1:length(phiseas), ")", sep = ""))
                
                if (length(thetaseasp) > 0) 
                  nm <- c(nm, paste("theta_p.", period, "(", 1:length(thetaseas), ")", sep = ""))
            } else {
                if (length(phi) > 0) 
                  nm <- c(nm, paste("phi(", 1:length(phi), ")", sep = ""))
                
                if (length(theta) > 0) 
                  nm <- c(nm, paste("theta(", 1:length(theta), ")", sep = ""))
                
                if (length(phiseas) > 0) 
                  nm <- c(nm, paste("phi.", period, "(", 1:length(phiseas), ")", sep = ""))
                
                if (length(thetaseas) > 0) 
                  nm <- c(nm, paste("theta.", period, "(", 1:length(thetaseas), ")", sep = ""))
            }
            
            if (length(dfrac) > 0) 
                nm <- c(nm, paste("d.f", sep = ""))
            if (length(H) > 0) 
                nm <- c(nm, paste("H", sep = ""))
            if (length(alpha) > 0) 
                nm <- c(nm, paste("alpha", sep = ""))
            
            if (length(dfs) > 0) 
                nm <- c(nm, paste("d.f.", period, sep = ""))
            if (length(Hs) > 0) 
                nm <- c(nm, paste("H.", period, sep = ""))
            if (length(alphas) > 0) 
                nm <- c(nm, paste("alpha.", period, sep = ""))
            
            if (itmean) {
              nm <- c(nm, "It. Fit. mean")
            } else  {
              if (is.logical(dmean) && !dmean) 
                nm <- c(nm, "zbar") else if (is.logical(dmean) && dmean) {
                  if(ans$strReg && !ans$intindex)
                    nm <- c(nm, "Intercept")
                  else if(!ans$intindex)
                    nm <- c(nm, "Fitted mean")
                  else
                    nm <- c(nm, ans$intname)
                }
              else nm <- c(nm, "Set mean")
            }
            if(ans$strReg){ 
              nm <- c(nm, names(xreglist$omega))
            }
            else {
              if (length(omega) > 0) {
                  for (j in 1:length(s)) nm <- c(nm, names(xreglist$omega[[j]]))
              }
              if (length(delta) > 0) {
                  for (j in 1:length(r)) nm <- c(nm, names(xreglist$delta[[j]]))
              }
              
            }
            
        }
        
        if (tpacf) 
            coeff <- c(phip, thetap, phiseasp, thetaseasp, dfrac, H, alpha, dfs, Hs, alphas, 
                       muHat, omega, delta)
        else
            coeff <- c(phi, theta, phiseas, thetaseas, dfrac, H, alpha, dfs, Hs, alphas, 
                       muHat, omega, delta)
        
        coeffs <- rbind(coeffs, coeff)
    }
    colnames(coeffs) <- nm
    rownames(coeffs) <- paste("Mode", 1:m)
    coeffs
}


AIC.arfima <- function(object, ..., k = 2) {
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    
    logl <- logLik(ans, ...)
    AIC(logl, k = k)
}

BIC.arfima <- function(object, ...) {
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    
    logl <- logLik(ans, ...)
    AIC(logl, k = log(length(object$z)))
}


vcov.arfima <- function(object, type = c("b", "o", "e"), cor = FALSE, digits = max(4, getOption("digits") - 
    3), tapprox = FALSE, summ = FALSE, ...) {
    if (!object$weeded) {
        warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
    }
    ans <- object
    
    if (type == c("b", "o", "e") && !summ) 
        cat("Returning both observed and theoretical covariance matrices: use type = \"o\" for just\n\t observed and \"e\" for just expected (from the expected Fisher information)\n")
    
    type <- type[1]
    
    xreg <- ans$xreg
    dmean <- ans$dmean
    itmean <- ans$itmean
    
    m <- length(ans$modes)
    
    ret <- vector("list", m)
    
    nm <- getnames(ans)
    
    if ((is.logical(ans$dmean) && !ans$dmean) || is.double(ans$dmean)) 
        nm <- nm[-length(nm)]
    period <- ans$period
    warnH <- warnX <- FALSE
    if (type %in% c("b", "o")) {
        for (i in 1:m) {
            ret[[i]]$observed <- solve(-ans$modes[[i]]$hess)
            
            if (cor) {
                temp <- ret[[i]]$observed
                
                for (j in 1:nrow(temp)) {
                  for (k in 1:ncol(temp)) {
                    temp[j, k] <- ret[[i]]$observed[j, k]/sqrt(ret[[i]]$observed[j, j] * 
                      ret[[i]]$observed[k, k])
                  }
                }
                ret[[i]]$observed <- temp
            }
            ret[[i]]$observed <- signif(ret[[i]]$observed, digits = digits)
            rownames(ret[[i]]$observed) <- colnames(ret[[i]]$observed) <- nm
        }
    }
    if (type %in% c("b", "e")) {
        if ((length(ans$modes[[1]]$H) > 0 || length(ans$modes[[1]]$Hs) > 0 || length(ans$modes[[1]]$alpha) > 
            0 || length(ans$modes[[1]]$alphas) > 0)) {
            if (!summ) 
                warning("The expected Fisher Information matrix has no known closed form for FGN or HD processes.\n\n\t\t\t\tReporting this matrix for the non-FGN/non-HD processes: however, results may be skewed")
            warnH <- TRUE
        }
        if (!is.null(xreg)) {
            if (!summ) 
                warning("The expected Fisher Information matrix does not include regression parameters.\n\n\t\t\t\tReporting this matrix for the ARFIMA processes: however, results may be skewed")
            warnX <- TRUE
        }
        ifdf <- length(ans$modes[[1]]$dfrac) > 0
        ifdfs <- length(ans$modes[[1]]$dfs) > 0
        for (i in 1:m) {
            cmode <- ans$modes[[i]]
            
            ret[[i]]$expected <- solve(iARFIMA(phi = cmode$phi, theta = cmode$theta, phiseas = cmode$phiseas, 
                thetaseas = cmode$thetaseas, dfrac = ifdf, dfs = ifdfs, period = period, 
                exact = !tapprox))/ans$n
            if (cor) {
                temp <- ret[[i]]$expected
                for (j in 1:nrow(temp)) {
                  for (k in 1:ncol(temp)) {
                    temp[j, k] <- ret[[i]]$expected[j, k]/sqrt(ret[[i]]$expected[j, j] * 
                      ret[[i]]$expected[k, k])
                  }
                }
                ret[[i]]$expected <- temp
            }
            ret[[i]]$expected <- signif(ret[[i]]$expected, digits = digits)
        }
    }
    names(ret) <- paste("Mode ", 1:m, sep = "")
    
    if (summ) {
        ret$warnH <- warnH
        ret$warnX <- warnX
    }
    
    ret
}


getnames <- function(ans) {
    
    if (class(ans) != "arfima") 
        stop("the getnames function only works for arfima objects")
    
    nm <- colnames(coef(weed(ans)))
    
    nm
} 
