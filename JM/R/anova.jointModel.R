anova.jointModel <-
function (object, object2, test = TRUE, 
        process = c("both", "Longitudinal", "Event"), L = NULL, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!missing(object2)) {
        if (!inherits(object2, "jointModel"))
            stop("Use only with 'jointModel' objects.\n")
        if (object$conv)
            warning("it seems that '", deparse(substitute(object)), 
                "' has not converged.\n")
        if (object2$conv)
            warning("it seems that '", deparse(substitute(object2)), 
                "' has not converged.")
        if (test && object$method != object2$method)
            stop("you compare joint models with different survival submodels.\n")
        L0 <- logLik(object)
        L1 <- logLik(object2)
        nb0 <- attr(L0, "df")
        nb1 <- attr(L1, "df")
        df <- nb1 - nb0
        if (test && df < 0)
            stop("'object' should be nested in 'object2'.\n")
        out <- list(nam0 = deparse(substitute(object)), L0 = L0, 
            aic0 = AIC(object), bic0 = AIC(object, k = log(attr(L0, "nobs"))), 
            nam1 = deparse(substitute(object2)), L1 = L1, aic1 = AIC(object2), 
            bic1 = AIC(object2, k = log(attr(L1, "nobs"))), df = df, test = test)
        if (test) {
            LRT <- - 2 * (L0 - L1)
            attributes(LRT) <- NULL
            if (LRT < 0)
                warning("either the two models are not nested or the model ", 
                    "represented by 'object2' fell on a local maxima.\n")
            out$LRT <- LRT
            out$p.value <- pchisq(LRT, df, lower.tail = FALSE)
        }
    } else {
        process <- match.arg(process)
        f <- function (t, thetas, V) {
            df <- length(t)
            L <- matrix(0, df, length(thetas))
            col.ind <- t
            L[cbind(seq_along(t), col.ind)] <- 1
            Ltheta <- c(L %*% thetas)
            LVtL <- L %*% tcrossprod(V, L)
            stat <- c(crossprod(Ltheta, solve(LVtL, Ltheta)))
            pval <- pchisq(stat, df, lower.tail = FALSE)
            list(Chisq = stat, df = df, "Pr(>|Chi|)" = pval)
        }
        VarCov <- vcov(object)
        if (process == "Longitudinal" || process == "both") {
            bts <- object$coefficients$betas
            indY <- seq(1, length(bts))
            terms <- object$assignY
            termsLabs <- all.vars(delete.response(object$termsYx))
            indMI <- lapply(termsLabs, 
                function (x) grep(x, names(terms), fixed = TRUE))
            indI <- lapply(grep(":", names(terms), fixed = TRUE), c)
            resY <- t(sapply(c(indMI, indI), function (tt) {
                xx <- unlist(terms[tt])
                f(xx, thetas = bts, V = VarCov[indY, indY])
            }))
            rownames(resY) <- c(termsLabs, names(terms)[unlist(indI)])
            resY <- as.data.frame(resY)    
        }
        if (process == "Event" || process == "both") {
            ga <- fixef(object, process = "Event")
            indT <- head(grep("T.", colnames(VarCov), fixed = TRUE), 
                length(ga))
            jj <- grep("Assoct[!^\\.s]", names(ga))
            ii <- setdiff(grep("Assoct", names(ga)), jj)
            termsLabs <- unique(c(names(object$assignT),
                all.vars(object$interFac$value), all.vars(object$interFac$slope)))
            AssoctLabs <- if (Asc.All <- length(object$coefficients$alpha) + 
                length(object$coefficients$Dalpha) == 1) {
                c("Assoct(?![.]s)", "Assoct([.]s)")
            } else {
                c("Assoct", "Assoct(?![.]s)", "Assoct([.]s)")
            }
            termsLabs2 <- switch(object$parameterization,
                "both" = if (Asc.All) 
                    c(termsLabs2, "Assoct", "Assoct.s") else 
                        c(termsLabs, "Assoct(all)", "Assoct", "Assoct.s"), 
                "value" = if (Asc.All) c(termsLabs, "Assoct") else
                    c(termsLabs, "Assoct(all)", "Assoct"),
                "slope" = if (Asc.All) c(termsLabs, "Assoct.s") else 
                    c(termsLabs, "Assoct(all)", "Assoct.s"))   
            indInt <- grep(":", termsLabs, fixed = TRUE)
            tt <- as.list(termsLabs)
            if (length(indInt))
                tt[indInt] <- strsplit(termsLabs[indInt], ":")
            indM <- lapply(tt, function (x) {
                if (length(x) > 1)
                    do.call(intersect, lapply(x, 
                        function (y) grep(y, names(ga), perl = TRUE)))
                else 
                    grep(x, names(ga), perl = TRUE)
            })
            indA <- lapply(AssoctLabs, 
                function (x) grep(x, names(ga), perl = TRUE))
            allEf <- c(indM, indA)
            allEf <- allEf[sapply(allEf, length) > 0]
            resT <- t(sapply(allEf, function (z) { 
                f(z, thetas = ga, V = VarCov[indT, indT])
            }))
            rownames(resT) <- termsLabs2
            resT <- as.data.frame(resT)
        }
        out <- if (is.null(L)) {
            switch(process, 
                "both" = list(aovTab.Y = resY, aovTab.T = resT),
                "Longitudinal" = list(aovTab.Y = resY),
                "Event" = list(aovTab.T = resT))
        } else {
            thetas <- switch(process, "both" = c(bts, ga),
                "Longitudinal" = bts, "Event" = ga)
            V <- switch(process, 
                "both" = VarCov[c(indY, indT), c(indY, indT)],
                "Longitudinal" = VarCov[indY, indY], 
                "Event" = VarCov[indT, indT])
            if (!is.numeric(L) || ncol(L) != length(thetas))
                stop("L matrix not of appropriate dimensions. ", 
                    "It should have ", length(thetas), " columns.\n")
            Ltheta <- c(L %*% thetas)
            LVtL <- L %*% tcrossprod(V, L)
            stat <- c(crossprod(Ltheta, solve(LVtL, Ltheta)))
            pval <- pchisq(stat, nrow(L), lower.tail = FALSE)
            res <- data.frame(Chisq = stat, df = nrow(L), 
                "Pr(>|Chi|)" = pval, check.names = FALSE, row.names = "L")
            out <- list(aovTab.L = res)
        }
    }
    class(out) <- "aov.jointModel"
    out
}
