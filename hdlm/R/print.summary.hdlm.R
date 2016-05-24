print.summary.hdlm <-
function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
      signif.stars= getOption("show.signif.stars"),...)
{
    cat("\nCall:\n", # S has ' ' instead of '\n'
    paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    resid <- x$residuals
    df <- x$rank
    rdf <- x$rank[1]
    if (rdf > 5L) {
        cat("Residuals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    zz <- zapsmall(quantile(resid), digits + 1)
    rq <- structure(zz, names = nam)
        
    print(rq, digits = digits, ...)
    } else if (rdf > 0L) {
        print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
        cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!\n")
    }

    z <- x

      if(z$level == 1 & z$pval.method != 'median') {
        index <- union(which(names(z$coefficients) == '(Intercept)'),
            intersect(which(z$p.value != 1),which(z$coefficients != 0)))
      } else if(z$level == 1 & z$pval.method == 'median') {
        index <- union(which(names(z$coefficients) == '(Intercept)'), which(z$coefficients != 0))
      } else if(z$level == 2) {
        index <- union(union(which(names(z$coefficients) == '(Intercept)'), which(z$coefficients != 0)), which(z$p.value != 1))
      } else {
        index <- 1:length(z$coefficients)
      } 
      z$lower.bound <- z$lower.bound[index]
      z$upper.bound <- z$upper.bound[index]
      z$p.value <- z$p.value[index]
      z$coefficients <- z$coefficients[index]


    cat("\nCoefficients:\n")
    if(is.null(z$p.value)) {
        MAT <- cbind(z$coefficients)
        colnames(MAT) <- c("Estimate")
        HDprintCoefmat(MAT, digits=digits, signif.stars=FALSE, na.print="NA", pval=FALSE)
    } else {
        MAT <- cbind(z$coefficients, z$lower.bound, z$upper.bound, z$p.value)
        colnames(MAT) <- c("Estimate", "Lower Bd", "Upper Bd", "P-value")
        HDprintCoefmat(MAT, digits=digits, signif.stars=signif.stars, na.print="NA")
    }

 
    cat("\nEstimated sigma: ",
    format(signif(x$sigma.hat, digits)), "", sep='')
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-squared:", formatC(x$r.squared, digits=digits))
        cat(",\tAdjusted R-squared:",formatC(x$adj.r.squared,digits=digits),
        "\nF-statistic:", formatC(x$fstatistic[1L], digits=digits),
        "on", x$fstatistic[2L], "and",
        x$fstatistic[3L], "DF,  p-value:",
        format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
                           x$fstatistic[3L], lower.tail = FALSE), digits=digits),
        "\n")
    }
    cat("\n")#- not in S
    invisible(x)
}

