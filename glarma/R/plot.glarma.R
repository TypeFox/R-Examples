plot.glarma <- function(x, which = c(1L,3L,5L,7L,8L,9L), fits = 1L:3L,
                        ask = prod(par("mfcol")) <
                              length(which) && dev.interactive(),
                        lwdObs = 1, lwdFixed = 1, lwdGLARMA = 1,
                        colObs = "black", colFixed = "blue", colGLARMA = "red",
                        ltyObs = 2, ltyFixed = 1, ltyGLARMA = 1,
                        pchObs = 1, legend = TRUE, residPlotType = "h",
                        bins = 10, line = TRUE, colLine = "red",
                        colHist = "royal blue", lwdLine = 2,
                        colPIT1 = "red", colPIT2 = "black",
                        ltyPIT1 = 1, ltyPIT2 = 2, typePIT = "l",
                        ltyQQ = 2, colQQ = "black",
                        titles, ...)
{
    ## ask = prod(par("mfcol")) < length(which) && dev.interactive()
    show <- rep(FALSE, 10)
    show[which] <- TRUE
    showFits <- rep(FALSE, 3)
    showFits[fits] <- TRUE
    ## construct list containing either NULLs or requested titles
    if (missing(titles)) {
        ## all titles will be default titles
        titles  <- vector("list", 10)
        titles[[5]] <- "Histogram of Uniform PIT"
        titles[[6]] <- "Q-Q Plot of Uniform PIT"
        titles[[7]] <- "Histogram of Randomized Residuals"
        titles[[8]] <- "Q-Q Plot of Randomized Residuals"
        titles[[9]] <- "ACF of Randomized Residuals"
        titles[[10]] <- "PACF of Randomized Residuals"
    } else {
        defaultTitles <- vector("list", 10)
        ## replace NULLs with requested titles
        defaultTitles[which] <- titles
        titles <- defaultTitles
    }

    if (x$type == "Poi" | x$type == "NegBin") {
        if (show[1L] & any(showFits == TRUE)) {
            obs.logic <- showFits[1]
            fixed.logic <- showFits[2]
            glarma.logic <- showFits[3]
            fits <- list(obs = x$y, fixed = exp(x$eta), glarma = x$mu)
            legendNames <- names(fits)
            titleNames <- c("Observed", "Fixed", "GLARMA")
            ltyAll <- c(ltyObs, ltyFixed, ltyGLARMA)
            lwdAll <- c(lwdObs, lwdFixed, lwdGLARMA)
            colAll <- c(colObs, colFixed, colGLARMA)
            yRange <- c(diff(range(fits$obs)), diff(range(fits$fixed)),
                        diff(range(fits$glarma)))
            yRange[!showFits] <- -Inf
            yLim <- which.max(yRange)
        }
        if (any(show[2L:4L] == TRUE)) {
            residuals <- x$residuals
            if (is.null(titles[[2]])) {
                titles[2] <- paste("ACF of", x$residType, "Residuals")
            }
            if (is.null(titles[[3]])) {
                titles[3] <- paste(x$residType, "Residuals")
            }
            if (is.null(titles[[4]])) {
                titles[4] <-
                    paste("Normal Q-Q Plot of ",x$residType, "Residuals")
            }
        }
        if (ask) {
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))
        }
        if (show[1L] & any(showFits == TRUE)) {
            dev.hold()
            if (is.null(titles[[1]])){
                main <- paste(titleNames[showFits], collapse = " vs ")
            } else {
                main <- titles[1]
            }
            ts.plot(fits[[yLim]], ylab = "Counts", xlab = "Time",
                    col = NA, main = main, ...)
            if (obs.logic)
                lines(fits$obs, lwd = lwdObs, lty = ltyObs, col = colObs)

            if (fixed.logic)
                lines(fits$fixed, lwd = lwdFixed, lty = ltyFixed,
                      col = colFixed)

            if (glarma.logic)
                lines(fits$glarma, lwd = lwdGLARMA, lty = ltyGLARMA,
                      col = colGLARMA)

            if (legend & any(showFits == TRUE)) {
                par(xpd = NA)
                mfrow <- par("mfrow")
                graph.param <-
                    legend("top", legend = legendNames[showFits],
                           lty = ltyAll[showFits],
                           ncol = 3,
                           cex = 0.7 - (mfrow[1] - 1)/10 - (mfrow[2] - 1)/10,
                           bty = "n", plot = FALSE)
                legend(graph.param$rect$left,
                       graph.param$rect$top + graph.param$rect$h,
                       legend = legendNames[showFits], col = colAll[showFits],
                       lwd = lwdAll[showFits], lty = ltyAll[showFits],
                       ncol = 3,
                       cex = 0.7 - (mfrow[1] - 1)/10 - (mfrow[2] - 1)/10,
                       bty = "n",
                       text.font = 4)
                par(xpd = FALSE)
            }
            dev.flush()
        }
        if (show[2L]) {
            dev.hold()
            acf(residuals, main = titles[2], ...)
            dev.flush()
        }
        if (show[3L]) {
            dev.hold()
            plot.ts(residuals, ylab = "Residuals", type = residPlotType,
                    main = titles[3], ...)
            dev.flush()
        }
        if (show[4L]) {
            dev.hold()
            qqnorm(residuals, ylab = "Residuals", main = titles[4], ...)
            abline(0, 1, lty = 2)
            dev.flush()
        }
        if (show[5L]) {
          dev.hold()
          histPIT(x, bins = bins, line = line, colLine = colLine,
                  colHist = colHist, lwdLine = lwdLine, main = titles[[5]], ...)
          dev.flush()
        }
        if (show[6L]) {
          dev.hold()
          qqPIT(x, bins = bins, col1 = colPIT1, col2 = colPIT2,
                lty1 = ltyPIT1, lty2 = ltyPIT2, type = typePIT,
                main = titles[[6]], ...)
          dev.flush()
        }
        if (any(show[7L:10L] == TRUE)) {
            rt <- normRandPIT(x)$rt
        }
        if (show[7L]) {
          dev.hold()
          hist(rt, breaks = bins, main = titles[[7]],
               col = colHist, xlab = expression(r[t]), ...)
          box()
          dev.flush()
        }
        if (show[8L]) {
          dev.hold()
          qqnorm(rt, main =  titles[[8]], ...)
          abline(0, 1, lty = ltyQQ, col = colQQ, ...)
          dev.flush()
        }
        if (show[9L]) {
          dev.hold()
          acf(rt,  main =  titles[[9]], ...)
          dev.flush()
        }
        if (show[10L]) {
          dev.hold()
          pacf(rt,  main =  titles[[10]], ...)
          dev.flush()
        }
    }

    if (x$type == "Bin") {
        if (show[1L] & any(showFits = TRUE)) {
            obs.logic <- showFits[1]
            fixed.logic <- showFits[2]
            glarma.logic <- showFits[3]
            observed <- x$y[, 1]/apply(x$y, 1, sum)
            fits <- list(fixed = 1/(1 + exp(-x$eta)),
                         glarma = 1/(1 + exp(-x$W)))
            legendNames <- c("obs", names(fits)[1], names(fits)[2])
            titleNames <- c("Observed", "Fixed", "GLARMA")
            pchAll <- c(pchObs, NA, NA)
            ltyAll <- c(ltyObs, ltyFixed, ltyGLARMA)
            lwdAll <- c(NA, lwdFixed, lwdGLARMA)
            colAll <- c(colObs, colFixed, colGLARMA)
        }
        if (any(show[2L:4L] == TRUE)) {
            residuals <- x$residuals
            residuals <- x$residuals
            if (is.null(titles[[2]])) {
                titles[2] <- paste("ACF of", x$residType, "Residuals")
            }
            if (is.null(titles[[3]])) {
                titles[3] <- paste(x$residType, "Residuals")
            }
            if (is.null(titles[[4]])) {
                titles[4] <-
                    paste("Normal Q-Q Plot of ",x$residType, "Residuals")
            }
        }
        if (ask) {
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))
        }
        if (show[1L] & any(showFits == TRUE)) {
            dev.hold()
            if (is.null(titles[[1]])){
                main <- paste(titleNames[showFits], collapse = " vs ")
            } else {
                main <- titles[1]
            }
            plot(1:length(observed), observed, ylab = "Counts", xlab = "Time",
                 col = NA, main = main, ...)
            if (obs.logic)
                points(observed, pch = pchObs, col = colObs)

            if (fixed.logic)
                lines(fits$fixed, lwd = lwdFixed, lty = ltyFixed,
                      col = colFixed)

            if (glarma.logic)
                lines(fits$glarma, lwd = lwdGLARMA, lty = ltyGLARMA,
                      col = colGLARMA)

            if (legend & any(showFits == TRUE)) {
                par(xpd = NA)
                mfrow <- par("mfrow")
                graph.param <-
                    legend("top", legend = legendNames[showFits],
                           pch = pchAll[showFits],
                           lty = ltyAll[showFits], ncol = 3,
                           cex = 0.7 - (mfrow[1] - 1)/10 - (mfrow[2] - 1)/10,
                           text.font = 4, plot = FALSE)

                legend(graph.param$rect$left,
                       graph.param$rect$top + graph.param$rect$h,
                       legend = legendNames[showFits], pch = pchAll[showFits],
                       col = colAll[showFits], lwd = lwdAll[showFits],
                       lty = ltyAll[showFits], ncol = 3,
                       cex = 0.7 - (mfrow[1] - 1)/10 - (mfrow[2] - 1)/10,
                       bty = "n", text.font = 4)
                par(xpd = FALSE)
            }
            dev.flush()
        }
        if (show[2L]) {
            dev.hold()
            acf(residuals, main = titles[2], ...)
            dev.flush()
        }
        if (show[3L]) {
            dev.hold()
            plot.ts(residuals, ylab = "Residuals", type = residPlotType,
                    main = titles[3], ...)
            dev.flush()
        }
        if (show[4L]) {
            dev.hold()
            qqnorm(residuals, ylab = "Residuals", main = titles[4], ...)
            abline(0, 1, lty = 2)
            dev.flush()
        }
        if (show[5L]) {
           dev.hold()
           histPIT(x, bins = bins, line = line, colLine = colLine,
                   colHist = colHist, lwdLine = lwdLine, main = titles[5], ...)
           dev.flush()
        }
        if (show[6L]) {
           dev.hold()
           qqPIT(x, bins = bins, col1 = colPIT1, col2 = colPIT2,
                 lty1 = ltyPIT1, lty2 = ltyPIT2, type = typePIT,
                 main = titles[6], ...)
           dev.flush()
        }
        if (any(show[7L:10L] == TRUE)) {
            rt <- normRandPIT(x)$rt
        }
        if (show[7L]) {
          dev.hold()
          hist(rt, breaks = bins, main = titles[[7]],
               col = colHist, xlab = expression(r[t]), ...)
          box()
          dev.flush()
        }
        if (show[8L]) {
          dev.hold()
          qqnorm(rt, main =  titles[[8]], ...)
          abline(0, 1, lty = ltyQQ, col = colQQ, ...)
          dev.flush()
        }
        if (show[9L]) {
          dev.hold()
          acf(rt,  main =  titles[[9]], ...)
          dev.flush()
        }
        if (show[10L]) {
          dev.hold()
          pacf(rt,  main =  titles[[10]], ...)
          dev.flush()
        }
    }
    invisible()
}
