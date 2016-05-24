
crosstab <- function(dep, indep, weight = NULL,
                     digits = list(expected = 1, prop = 3, percent = 1, others = 3),
                     max.width = NA, expected = FALSE, prop.r = FALSE,
                     prop.c = FALSE, prop.t = FALSE, prop.chisq = FALSE,
                     chisq = FALSE, fisher = FALSE, mcnemar = FALSE, resid = FALSE,
                     sresid = FALSE, asresid = FALSE, missing.include = FALSE,
                     drop.levels = TRUE, format = "SPSS", cell.layout = TRUE,
                     row.labels = !cell.layout,
                     percent = (format == "SPSS" && !row.labels),
                     total.r, total.c,
                     dnn = NULL, xlab = NULL, ylab = NULL, main = "",
                     user.missing.dep, user.missing.indep,
                     plot = getOption("descr.plot"), ...)
{
    if(missing(dep))
        stop("The argument 'dep' (dependent variable) is missing.")
    if(missing(indep))
        stop("The 'indep' (independent variable) is missing. Please, consider using either CrossTable() or freq().")

    if(is.null(dnn))
        dnn <- c(deparse(substitute(dep)), deparse(substitute(indep)))

    if(!missing(user.missing.indep)){
        user.missing.indep <- paste("^", user.missing.indep, "$", sep = "")
        ilevels <- levels(indep)
        for(lev in user.missing.indep){
            if(length(grep(lev, ilevels))){
                idx <- grep(lev, as.character(indep))
                if(length(idx))
                    indep[idx] <- NA
            }
        }
        indep <- factor(indep)
    }
    if(!missing(user.missing.dep)){
        user.missing.dep <- paste("^", user.missing.dep, "$", sep = "")
        dlevels <- levels(dep)
        for(lev in user.missing.dep){
            if(length(grep(lev, dlevels))){
                idx <- grep(lev, as.character(dep))
                if(length(idx))
                    dep[idx] <- NA
            }
        }
        dep <- factor(dep)
    }
    if(missing.include){
        dep <- no.drop.levels(dep)
        indep <- no.drop.levels(indep)
    }
    if(drop.levels){
        dep <- factor(dep)
        indep <- factor(indep)
    }
    if (is.null(weight))
        tab <- table(dep, indep)
    else
        tab <- round(xtabs(weight ~ dep + indep))
    names(dimnames(tab)) <- dnn

    if(!missing(total.r)){
        if(!is.logical(total.r))
            stop(gettext("total.r must be logical", domain = "R-descr"))
        if(missing(total.c))
            total.c <- total.r
    }
    if(!missing(total.c)){
        if(!is.logical(total.c))
            stop(gettext("total.c must be logical", domain = "R-descr"))
        if(missing(total.r))
            total.r <- total.c
    }
    if(missing(total.r) & missing(total.c))
        total.r <- total.c <- TRUE

    crosstb <- CrossTable(tab, digits = digits, max.width = max.width,
                          expected = expected, prop.r = prop.r,
                          prop.c = prop.c, prop.t = prop.t,
                          prop.chisq = prop.chisq, chisq = chisq,
                          fisher = fisher, mcnemar = mcnemar, resid = resid,
                          sresid = sresid, asresid = asresid,
                          missing.include = missing.include,
                          drop.levels = drop.levels, format = format, dnn = dnn,
                          cell.layout = cell.layout, row.labels = row.labels,
                          percent = percent, total.r = total.r,
                          total.c = total.c, xlab = xlab, ylab = ylab)

    if(plot == TRUE)
        plot.CrossTable(crosstb, ...)

    crosstb
}


plot.CrossTable <- function(x, xlab, ylab, main = "", col, inv.x = FALSE, inv.y = FALSE, ...)
{
    tabforplot <- t(x$tab)
    if(missing(xlab)){
        lab <- attr(x, "xlab")
        if(is.null(lab))
            xlab <- x$ColData
        else
            xlab <- lab
    }
    if(missing(ylab)){
        lab <- attr(x, "ylab")
        if(is.null(lab))
            ylab <- x$RowData
        else
            ylab <- lab
    }
    nxlev <- dim(tabforplot)[1]
    nylev <- dim(tabforplot)[2]
    if(missing(col)){
        col.min <- 0.9 - 0.25 * (nylev - 1)
        if(col.min < 0.3)
            col.min  <- 0.3
        col <- gray.colors(nylev, 0.9, col.min)
    }
    if(inv.x)
        tabforplot <- tabforplot[nxlev:1, ]
    if(inv.y)
        tabforplot <- tabforplot[, nylev:1]
    class(tabforplot) <- "table"
    if(length(grep("^color$", names(list(...)))) == 0)
        mosaicplot(tabforplot, main = main, xlab = xlab, ylab = ylab, col = col, ...)
    else
        mosaicplot(tabforplot, main = main, xlab = xlab, ylab = ylab, ...)
}

