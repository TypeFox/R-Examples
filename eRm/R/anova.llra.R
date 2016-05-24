anova.llra <- function(object, ...) UseMethod("anova.llra")

anova.llra <- function(object,...)
  {
    objets <- list(object, ...)
    isllra <- unlist(lapply(objets, function(x) "llra" %in% class(x)))

    ## checks
    if (!all(isllra)) {
        objets <- objets[isllra]
        warning("non-LLRA-model(s) removed")
    }
    dimdata <- dim(objets[[1L]]$X)
    samedata <- unlist(lapply(objets, function(x) all(dim(x$X)==dimdata)))
    if (!all(samedata))
        stop("models were not all fitted to the same size of dataset")

    nmodels <- length(objets)
    logliks <- as.numeric(lapply(objets, function(x) x$loglik))
    npars <- as.numeric(lapply(objets, function(x) x$npar))
    nparsS <- npars[order(npars)]
    logliksS <- logliks[order(npars)]
    lrstat <-c(NA,2*abs(diff(logliksS)))
    dfs <- c(NA,abs(diff(nparsS)))
    ps <- 1-pchisq(lrstat,dfs)
    tbl <- data.frame(nparsS, logliksS, dfs, lrstat, ps)
    dimnames(tbl) <- list(1L:nmodels, c("Npar", "logLik", "df", "-2LR","Asymp.p-Value"))
    title <- "Analysis of Deviance Table\n"
    structure(tbl, heading = title, class = c("anova",
        "data.frame"))
}
