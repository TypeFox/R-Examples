RFcluster <-
    function (formula, id, data, nfold = 15,
              ntree=500, progress=TRUE, printit = TRUE, seed = 29)
{
    m <- match.call(expand.dots = FALSE)
    idnam <- deparse(m[["id"]])
    names(m)[2] <- "formula"
    mm <- match(c("formula", "data", "na.action"), names(m), 0L)
    m <- m[c(1L, mm)]
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1L]] <- as.name("model.frame")
    m[[2]] <- as.formula(paste(deparse(m[[2]]),"-",idnam,sep=""))
    mf <- eval(m, parent.frame())
    Terms <- attr(mf, "terms")
    attr(Terms, "intercept") <- 0
    y <- model.response(mf)
    ylev <- levels(y)
    yfac <- factor(ylev, levels=ylev)
    mf <- model.frame(terms(reformulate(attributes(Terms)$term.labels)),
                      data.frame(mf))
    id <- data[,idnam]
    clusts <- levels(id)
    ynam <- all.names(formula)[1]
    sampfreqs <- prop.table(table(y))
    ncat <- length(clusts)
    mat <- matrix(0, length(id), ntree)
        for (i in 1:ntree) {
        if (progress & i%%10 == 0)
            cat(i, "")
        samp <- sample(clusts, replace = TRUE)
        inbag <- id %in% samp
        outbag <- !inbag
        pred <- randomForest(mf[inbag, ], y[inbag], xtest = mf[outbag,
                                                    ], ntree = 1,
                             ytest = y[outbag])$test$predicted
        mat[outbag, i] <- unclass(pred)
    }
    if (progress)
        cat("\n")
    mat[mat == 0] <- NA
    categ <- apply(mat, 1, function(x) {
        tab <- table(ylev[x])
        names(tab)[which.max(tab)]
    })
    categ <- factor(categ, levels=ylev)
    tab <- table(y, categ)
    accmat <- t(apply(tab, 1, function(x) x/sum(x)))
    acc <- sum(sampfreqs*diag(accmat))
    if (printit)
         cat("OOB accuracy =", round(acc, 2), "\\n")
    invisible(list(class=categ, OOBaccuracy=acc, confusion = accmat))
}
