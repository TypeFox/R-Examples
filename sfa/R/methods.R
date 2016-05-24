coef.sfa <- function(object, ...)
    return(object$coef)

print.sfa <- function(x, ...) {
    cat("\n")
    cat("Stochastic frontier analysis model\n")
    cat("\n")
    cat("Coefficients:\n")
    print(coef(x))
    cat("\n")
    invisible(x)
}

predict.sfa <- function(object, newdata = NULL, intercept = NULL, ...)  {
    cf <- coef(object)[1:(length(coef(object))-2)]
    if(object$fun == "tnormal" & is.null(object$par_mu)) cf <- coef(object)[1:(length(coef(object))-3)]
    if(is.null(intercept) && names(cf)[1] != "Intercept") intercept <- FALSE
    if(is.null(intercept) && names(cf)[1] == "Intercept") intercept <- TRUE
    if(names(cf)[1] != "Intercept" && intercept){
        warning("Intercept coef missing, setting intercept to FALSE",
                call. = FALSE, immediate. = TRUE)
        intercept <- FALSE
    }
    if(is.null(newdata)) {
        X <- object$X
    } else {
        X <- newdata
        if (intercept) {
            cbind(1, X)
        }
    }
    X%*%cf
}

fitted.sfa <- function(object, ...){
    predict(object, ...)
}

logLik.sfa <- function(object, ...)
    return(object$logLik)
    
residuals.sfa <- function(object, ...) {
    return(object$y - fitted(object))
}

summary.sfa <- function(object,...){
    coef <- coef(object)
    var_beta <- diag(solve(object$hess))
   	tvalue <- coef/sqrt(var_beta)
    coef.table <- cbind(coef, sqrt(var_beta), tvalue)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error","t value"))
    row.names(coef.table) <- names(object$coef)
    l_sfa <- object$logLik
    l_ols <- logLik(object$ols)
    L <- round(2*(l_sfa-l_ols)[1], digits = 3)
    names(L) <- "LR-Test"
    chisq_df <- attributes(logLik(object$ols))$df
    names(chisq_df) <- "df"
    p <- pchisq(L, chisq_df)
    names(p) <- "p-value"
    test.table <- rbind(l_sfa, l_ols)
    dimnames(test.table) <- list(NULL, "value Log-Lik")
    row.names(test.table) <- c("sfa", "ols")
    wert.table <- list(L=L, p=p, chisq_df=chisq_df)
    mean.eff <- mean(eff(object))
    names(mean.eff) <- "mean efficiency"
    ret <- list(coef.table=coef.table, test.table = test.table, p = p, chisq_df = chisq_df, LR = L, mean.eff = mean.eff)
### class
    class(ret) <- "summary.sfa"
    return(ret)
}

print.summary.sfa <- function(x, ...) {
    cat("Stochastic frontier analysis model\n")
    cat("\n")
    print(x$coef.table)
    cat("\n")
    cat("LR-test: sigmau2 = 0 (inefficiency has no influence to the model)")
    cat("\n")
    cat("H0: sigmau2 = 0 (beta_sfa = beta_ols)")
    cat("\n")
    cat("\n")
    print(x$test.table)
    cat("\n")
    cat(paste("value LR-Test:", x$LR))
    cat(paste(" on", x$chisq_df, "degrees of freedom"))
    cat(paste(" p-value", round(x$p, digits = 5)))
    cat("\n")
    cat("\n")
    print(x$mean.eff)
}

# nach Jodrow et al (1982)    
u.sfa <- function(object, ...) {
    sigmau2 <- object$sigmau2
    sigmav2 <- object$sigmav2
    fun <- object$fun
    sc <- object$sc
    if (fun == "hnormal") {
    mu_i <- -sc * residuals(object) * sigmau2 / (sigmau2 + sigmav2)
    sigma_i <- sqrt(sigmau2 * sigmav2 / (sigmau2 + sigmav2))
    }
    if (fun == "exp") {
    mu_i <- -sc * residuals(object) - sigmav2 / sqrt(sigmau2)
    sigma_i <- sqrt(sigmav2)
    }
    if (fun == "tnormal") {
    mu <- object$mu
    mu_i <- (-sc * residuals(object) * sigmau2 + mu * sigmav2) / (sigmau2 + sigmav2)
    sigma_i <- (sqrt(sigmau2) * sqrt(sigmav2)) / (sqrt(sigmau2 + sigmav2))
    }
    u <- mu_i + sigma_i*(dnorm(-mu_i/sigma_i)/pnorm(mu_i/sigma_i))
    return(u)
}

te.eff.sfa <- function(object, ...) {
    sigmau2 <- object$sigmau2
    sigmav2 <- object$sigmav2
    fun <- object$fun
    sc <- object$sc
    if (fun == "hnormal") {
    mu_i <- -sc * residuals(object) * sigmau2 / (sigmau2 + sigmav2)
    sigma_i <- sqrt(sigmau2) * sqrt(sigmav2) / sqrt(sigmau2 + sigmav2)
    }
    if (fun == "exp") {
    mu_i <- -sc * residuals(object) - sigmav2 / sqrt(sigmau2)
    sigma_i <- sqrt(sigmav2)
    }
    if (fun == "tnormal") {
    mu <- object$mu
    mu_i <- (-sc * residuals(object) * sigmau2 + mu * sigmav2) / (sigmau2 + sigmav2)
    sigma_i <- (sqrt(sigmau2) * sqrt(sigmav2)) / (sqrt(sigmau2 + sigmav2))
    }
    eff <- (1 - pnorm(sc * sigma_i - mu_i / sigma_i)) / (1 - pnorm(-mu_i / sigma_i))* exp(-sc * mu_i + 0.5 * sigma_i^2)
    return(eff)
    }

eff <- function(object, ...) {
    UseMethod("eff")
    }
    
# nach Coelli UNE:
eff.sfa <- function(object, ...) {
    sc <- object$sc
    eff <- (predict(object) - sc * u.sfa(object)) / predict(object)
    if (sc == -1) eff <- 1/eff
    return(eff)
    }
