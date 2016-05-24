HDprintCoefmat <-
function (x, digits = max(3, getOption("digits") - 2),
                            signif.stars = getOption("show.signif.stars"), 
                            signif.legend = signif.stars,
                            dig.tst = max(1, min(5, digits - 1)),
                            cs.ind = 1:nrow(x),
                            zap.ind = integer(0), 
                            P.values = NULL,
                            has.Pvalue = nc >= 4 && substr(colnames(x)[nc], 1, 3) == "Pr(",
                            eps.Pvalue = .Machine$double.eps,
                            na.print = "NA", pval = TRUE, ...) 
{
    if (is.null(d <- dim(x)) || length(d) != 2L) 
        stop("'x' must be coefficient matrix/data frame")
    nc <- d[2L]
    xm <- data.matrix(x)

    Cf <- array("", dim = d, dimnames = dimnames(xm))
    if(pval == TRUE) {
        Cf[,1:(nc-1)] <- format(x[,1:(nc-1)], digits=dig.tst)
    
            pv <- as.vector(xm[, nc])
            Cf[, nc] <- format.pval(pv, digits = dig.tst, 
                eps = eps.Pvalue)
            signif.stars <- signif.stars && any(pv < 0.1)
            
                Signif <- symnum(pv, corr = FALSE, na = FALSE, 
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  symbols = c("***", "**", "*", ".", " "))
                Cf <- cbind(Cf, format(Signif))
    } else {
        Cf[,1:(nc)] <- format(x[,1:(nc)], digits=dig.tst)
    }

    print.default(Cf, quote = FALSE, right = TRUE, na.print = na.print, 
        ...)

    if(pval == TRUE)  cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
    invisible(x)
}

