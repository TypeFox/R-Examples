CNVtest <-
function (x, type = "Wald")
{
    nCov <- attr(x, "nCov")
    F <- qr.solve(-x$hessian)
    type.test <- charmatch(type, c("Wald", "LRT"))
    cc <- NCOL(x$coefficients)
    if (type.test == 1) {
        if (attr(x, "model") == 1) {
            K <- diag(1, cc)[-cc, ] - diag(1, cc)[-1, ]
            beta <- x$coefficients[1, ]
            Var <- F[1:cc, 1:cc]
            stat <- as.double(t(K %*% beta) %*% qr.solve(K %*%
                Var %*% t(K)) %*% (K %*% beta))
            df <- nrow(K)
        }
        else {
            beta <- x$coefficients[2, 1]
            se <- sqrt(F[2, 2])
            stat <- (beta/se)^2
            df <- 1
        }
    }
    else {
        formula <- x$formula
        formula.null <- eval(parse(text = paste("update(formula,.~.-", x$CNVname, ")", sep = "")))
        family <- attr(x, "family")
        if (family != "weibull")
          model.null <- glm(formula.null, data = x$data, family = family)
        else
          model.null <- survreg(formula.null, data = x$data)          
        stat <- 2 * (logLik(x)[1] - logLik(model.null)[1])
        df <- if (attr(x, "model") == 1)
            cc - 1
        else 1
    }
    pvalue <- pchisq(stat, df, lower.tail = FALSE)
    out <- list(type = type.test, stat = stat, df = df, pvalue = pvalue)
    class(out) <- "CNVtest"
    out
}

