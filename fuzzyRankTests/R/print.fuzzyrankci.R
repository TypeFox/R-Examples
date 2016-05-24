print.fuzzyrankci <- function(x, digits = 4, ...)
{
    cat("\n")
    writeLines(strwrap(x$method, prefix = "\t"))
    cat("\n")
    cat("data: ", x$data.name, "\n")
    cat(format(100 * x$conf.level), "percent confidence interval:\n")
    cat("\n")

    if (x$conf.level == 0.0 || x$conf.level == 1.0) {
        writeLines(strwrap(paste("Fuzzy confidence interval has membership function that is", format(x$conf.level), "everywhere\n")))
    } else {

        xx <- x$knots
        vv <- x$knot.values
        ii <- x$interval.values
        if (is.finite(min(xx))) {
            xx <- c(-Inf, xx)
            vv <- c(NA, vv)
            ii <- c(0, ii)
        }
        if (is.finite(max(xx))) {
            xx <- c(xx, Inf)
            vv <- c(vv, NA)
            ii <- c(ii, 0)
        }

        foom <- ii[- length(ii)] + diff(ii) / 2
        barm <- vv[! is.na(vv)]
        ##### figure out whether regular type ####
        if (length(foom) == 4 && length(barm) == 4 &&
            all(abs(foom - barm) < x$tol) && x$interval.values[2] == 1) {
            writeLines(strwrap("Randomized confidence interval is mixture of two intervals\n"))
            cat("\n")
            foobaz <- matrix(c(x$interval.values[1], 1 - x$interval.values[1],
                x$knots[c(1, 2, 4, 3)]), nrow = 2)
            dimnames(foobaz) <- list(rep("", 2),
                c("probability", "lower end", "upper end"))
            print.default(foobaz, digits = digits, quote = FALSE, right = TRUE)
            cat("\n")
            writeLines(strwrap(paste("Corresponding fuzzy confidence interval is one on the narrower interval,", format(x$interval.values[1], digits = digits), "elsewhere on the wider interval, and zero outside the wider interval, with values at jumps that are the average of the left and right limits")))
        } else {
            writeLines(strwrap(paste("Fuzzy confidence interval:")))
            cat("\n")
            fred <- paste("{", format(xx[1], digits = digits), "}", sep = "")
            sally <- vv[1]
            for (i in 2:length(xx)) {
                fred <- c(fred,
                    paste("(", format(xx[i - 1], digits = digits), ", ",
                        format(xx[i], digits = digits), ")", sep = ""))
                sally <- c(sally, ii[i - 1])
                fred <- c(fred,
                    paste("{", format(xx[i], digits = digits), "}", sep = ""))
                sally <- c(sally, vv[i])
            }
            fred <- fred[! is.na(sally)]
            sally <- sally[! is.na(sally)]
            sally <- format(sally, digits = digits)
            foobaz <- cbind(fred, sally)
            dimnames(foobaz) <- list(rep("", nrow(foobaz)), c("set", "value"))
            print.default(foobaz, digits = digits, quote = FALSE, right = TRUE)
            cat("\n")
        }
    }

    cat("\n")
    invisible(x)
}
