cv.lm <-
function (data = DAAG::houseprices, form.lm = formula(sale.price ~ area),
              m = 3, dots = FALSE, seed = 29, plotit = c("Observed","Residual"),
              main="Small symbols show cross-validation predicted values",
              legend.pos="topleft", printit = TRUE)
{
    gphtype <- ""
    if(is.logical(plotit)){
        if(plotit)gphtype <- "Observed"
        } else
    if(is.character(plotit)){
             if(!(plotit[1]%in%c("Observed","Residual","")))
             stop (paste("Illegal argument plotit =", plotit[1]))
             gphtype <- plotit[1]
             if(plotit[1]%in%c("Observed","Residual"))plotit <- TRUE
         } else stop("Argument plotit must be logical or character")
    if (class(form.lm) == "formula") form <- form.lm else
    if (class(form.lm) %in% c("call","lm")) form <- formula(form.lm) else
    stop("form.lm must be formula or call or lm object")
    formtxt <- deparse(form)
    mf <- model.frame(form, data=data)
    ynam <- attr(mf, "names")[attr(attr(mf, "terms"),"response")]
    data.lm <- lm(mf)
    tm <- terms(mf)
    xcolumns <- labels(tm)
    n <- nrow(data)
    data[, ynam] <- model.response(mf)
    data[, "Predicted"] <- predict(data.lm)
    data[, "cvpred"] <- numeric(n)
    yval <- mf[, ynam]
    if(gphtype=="Residual") yval <- yval-data[, "Predicted"]
    if (!is.null(seed)) set.seed(seed)
    n <- dim(data)[1]
    rand <- sample(n)%%m + 1
    foldnum <- sort(unique(rand))
    for (i in foldnum) {
        rows.in <- rand != i
        rows.out <- rand == i
        subs.lm <- lm(form, data=data[rows.in, ])
        data[rows.out, "cvpred"] <- predict(subs.lm, newdata=data[rows.out,])
    }
    if (length(xcolumns)==1){
        stline <- TRUE
        xnam <- xcolumns
    } else {
        stline <- FALSE
        xnam <- "Predicted"
    }
    if (printit) {
        options(digits = 3)
        print(anova(data.lm))
        cat("\n")
    }
    if (plotit) {
        oldpar <- par(mar = par()$mar - c(1, 0, 2, 0))
        on.exit(par(oldpar))
        coltypes <- palette()[c(2, 3, 6, 1, 4:5, 7)]
        if (m > 7)
            coltypes <- c(coltypes, rainbow(m - 7))
        ltypes <- 1:m
        ptypes <- 2:(m + 1)
        par(lwd = 2)
        if(stline) xlab <- xnam else {
            xlab <- "Predicted (fit to all data)"
            cat("\n")
            warning(paste("\n\n As there is >1 explanatory variable, cross-validation\n",
                          "predicted values for a fold are not a linear function\n",
                          "of corresponding overall predicted values.  Lines that\n",
                          "are shown for the different folds are approximate\n"))
        }
        ylab <- ynam
        if(gphtype=="Residual")ylab <- paste(ynam, " (offset from predicted using all data)")
        plot(as.formula(paste("yval ~", xnam)), data=data, ylab = ylab, type = "p",
             pch=ptypes[rand], col=coltypes[rand], cex=1.25,
             xlab=xlab)
        title(main=main, cex=1.05)
        if(dots){
            with(data,
        points(as.formula(paste("yval ~", xnam)), data=data, type = "p",
             pch=16, col=coltypes[rand], cex=1)
                 )
    }
    }
    if(printit|plotit){
        sumss <- 0
        sumdf <- 0
        for (i in foldnum) {
            rows.in <- rand != i
            rows.out <- rand == i
            n.out <- sum(rows.out)
            resid <- data[rows.out, ynam] - data[rows.out, "cvpred"]
            ss <- sum(resid^2)
            sumss <- sumss+ss
            if (printit) {
                fold_data <- t(cbind(data[rows.out, c(xnam, "cvpred", ynam)], resid))
                rownames(fold_data) = c(xnam, "cvpred", ynam, "CV residual")
                cat("\nfold", i, "\n")
                cat("Observations in test set:", n.out, "\n")
                print(fold_data, collab = rep("", n.out))
                cat("\nSum of squares =", round(ss, 2), "   Mean square =",
                    round(ss/n.out, 2), "   n =", n.out, "\n")
            }
            if (plotit) {
                xval <- data[rows.out, xnam]
                nminmax <- c(which.min(xval), which.max(xval))
                cvpred <- data[rows.out, "cvpred"]
                if(gphtype=="Residual")cvpred <- cvpred-data[rows.out, "Predicted"]
                points(xval, cvpred, col = coltypes[i], pch = ptypes[i],
                       cex = 0.75, lwd=1)
                n1 <- which.min(xval)
                n2 <- which.max(xval)
                fold.lm <- lm(cvpred ~ xval)
                fold.b <- coef(fold.lm)
                lines(xval[c(n1,n2)], fold.b[1]+fold.b[2]*xval[c(n1,n2)],
                      col = coltypes[i], lty = ltypes[i])
                topleft <- par()$usr[c(1, 4)]
                par(lwd = 1, col=1)
                legend(x=legend.pos, legend = paste("Fold",
                                               1:m), pch = ptypes,
                       lty = ltypes, col = coltypes,
                       cex = 0.75)
            }
        }
    }
    sumdf <- sum(!is.na(data[, "Predicted"]))
    if(printit){
        cat("\nOverall", "(Sum over all", n.out, "folds)", "\n")
        print(c(`ms` = sumss/sumdf))
    }
    attr(data, "ms") <- sumss/sumdf
    attr(data, "df") <- sumdf
    invisible(data)
}
