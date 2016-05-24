lrtest.splm <- function(x, y, ...) {
    ## correspondence with component names
    compnames <- c("RE", "AR(1)", "SEM", "SAR")
    compnames.long <- c("random effects", "AR(1) errors",
                        "spatial errors", "spatial lag")
    names(compnames) <- c("phi", "rho", "lambda", "psi")

    ## first: check same betas!
    if(!identical(names(x$coef), names(y$coef))) {
        stop("Models estimated on different regressors")
    }
    ## ...and Nobs
    if(!identical(length(x$residuals), length(y$residuals))) {
        stop("Models estimated on different number of obs.")
    }

    ## build error and ev. SAR components vector
    ecompsx <- c(names(x$errcomp), names(x$arcoef))
    ecompsy <- c(names(y$errcomp), names(y$arcoef))

    ## check that models be different
    if(identical(ecompsx, ecompsy)) stop("The model is the same")

    ## which model is bigger?
    if(length(ecompsx)>=length(ecompsy)) {
        m1 <- x
        m0 <- y
        ecomps1 <- ecompsx
        ecomps0 <- ecompsy
    } else {
        m1 <- y
        m0 <- x
        ecomps1 <- ecompsy
        ecomps0 <- ecompsx
    }

    ## check if nested
    if(!all(ecompsx %in%ecompsy) & !all(ecompsy %in% ecompsx)) {
        stop("Models are not nested")
    }

    ll1 <- m1$logLik
    ll0 <- m0$logLik

    testedparms <- ecomps1[-which(ecomps1 %in% ecomps0)]

    LRstat <- 2*(ll1-ll0)
    df <- abs(length(ecompsx) - length(ecompsy))
    names(df) <- "df"
    pLR <- pchisq(LRstat, df = df, lower.tail=FALSE)

    names(LRstat) <- "LR test"
    RVAL <- list(statistic = LRstat, parameter = df,
                 method = paste("Likelihood ratio for exclusion of",
                 paste(compnames[testedparms], collapse = ", "), "\n\n",
                 "from panel model with ",
                 paste(compnames[ecomps1], collapse = ", ")),
                 p.value = pLR, data.name = "dname here (see bgtest)")
    class(RVAL) <- "htest"
    return(RVAL)
}
