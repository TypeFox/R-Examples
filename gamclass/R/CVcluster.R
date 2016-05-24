CVcluster <-
    function (formula, id, data, na.action=na.omit, nfold = 15, FUN = lda,
              predictFUN=function(x, newdata, ...)predict(x, newdata, ...)$class,
              printit = TRUE, cvparts = NULL, seed = 29)
{
    mf <- match.call(expand.dots = FALSE)
    idnam <- deparse(mf[["id"]])
    mm <- match(c("formula", "data", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, mm)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    formula <- mf[[2]]

    mf <- eval(mf, parent.frame())
    if(idnam%in%names(mf)){
        formtxt <- deparse(formula)
        if(grep('.', formtxt)){
            formtxt <- sub('.', paste("(.-", idnam, ")", sep=""),
                           formtxt, fixed=TRUE)
            formula <- as.formula(formtxt)
        }
        id <- mf[,idnam]
    }
    idval <- unique(id)
    if (is.null(cvparts)) {
        set.seed(seed)
        if (nfold == length(idval))
            repl <- FALSE
        else repl <- TRUE
        cvparts <- sample(1:nfold, length(idval), replace = repl)
    }
    y <- hat <- eval(formula[[2]], envir = as.data.frame(data))
    for (i in cvparts) {
        testclust <- idval[cvparts == i]
        testrows <- id %in% testclust
        trainrows <- !testrows
        model <- FUN(formula, data = data[trainrows, ])
        hat[testrows] <- predictFUN(model, newdata = data[testrows,
                                           ])
    }
    tab <- table(y, hat)
    accmat <- t(apply(tab, 1, function(x) x/sum(x)))
    acc <- sum(tab[row(tab) == col(tab)])/sum(tab)
    if (printit)
        cat("CV accuracy =", round(acc, 2), "\n")
    invisible(list(class = hat, CVaccuracy = acc, confusion = accmat))
}
