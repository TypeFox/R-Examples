diagBinom <- function(mod, type = 1:2,
                      ask = prod(par("mfcol")) < length(type) && dev.interactive()) {
    stopifnot(mod$family$family == "binomial")
    
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    fittedValues <- fitted(mod)
    y <- mod$y
    
    if (1 %in% type) {
        decFittedValues <- cut(fittedValues,
                               breaks = unique(quantile(fittedValues,
                                                        seq(0, 1, .1))),
                               include.lowest = TRUE)
        observedValues <- tapply(y, decFittedValues, mean)
        meanFittedValues <- tapply(fittedValues, decFittedValues, mean)
        plot(observedValues, meanFittedValues,
             xlim = 0:1, ylim = 0:1,
             main = "Observed probability\nvs. Fitted",
             xlab = "Observed", ylab = "Fitted")
        abline(a = 0, b = 1, lty = 2)
    }
    
    if (2 %in% type) {
        if (length(table(y)) > 2) {
            warning("Logistic regression estimated from aggregated data, ROC curve cannot be estimated.")
        } else {
            y2 <- y == max(y)
            totalPos <- sum(y2)
            totalNeg <- length(y2) - totalPos
            
            res <- sapply(sort(fittedValues), function(x) {
                se <- sum(y2[fittedValues > x]) / totalPos
                sp <- sum(! y2[fittedValues <= x]) / totalNeg
                c(se = se, sp = sp)
            })
            
            plot(c(0, res["sp", ], 1), c(1, res["se", ], 0), type = "l",
                 xlim = c(1, 0), ylim = c(0, 1),
                 xlab = "Specificity", ylab = "Sensitivity",
                 main = "ROC curve")
            abline(a = 1, b = -1, lty = 2)
        }
            
    }
}
