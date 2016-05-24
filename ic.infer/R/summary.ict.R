summary.ict <- function (object, brief = TRUE, digits = max(3, getOption("digits") - 
    3), scientific = FALSE, tol = sqrt(.Machine$double.eps), 
    ...) 
{
    if (!"ict" %in% class(object)) 
        stop("summary.ict works on class ict only.")
    cat("Order-related hypothesis test:", "\n")
    namen <- names(object$b.restr)
    if (is.null(namen)) 
        namen <- paste("m", 1:length(object$b.restr), sep = "")
    hilf <- object$ui
    colnames(hilf) <- namen[object$restr.index]
    aus.rest <- cbind(hilf, rep("%*%colnames", nrow(object$ui)), 
        c(rep(" == ", object$meq), rep(" >= ", nrow(object$ui) - 
            object$meq)), object$ci)
    rownames(aus.rest) <- paste(1:nrow(aus.rest), ":", sep = "")
    colnames(aus.rest)[(ncol(hilf) + 1):ncol(aus.rest)] <- rep(" ", 
        ncol(aus.rest) - ncol(hilf))
    aus.rest <- cbind(rep(" ", nrow(aus.rest)), aus.rest)
    aus.rest[object$iact, 1] <- "A"
    aus.test <- c(format(object$T, digits = digits), format(if (object$p.value < 
        1e-04) "<0.0001" else round(object$p.value, 4), digits = 4))
    names(aus.test) <- c("Test statistic", "p-value")
    if (object$TP == 1) {
        cat("\nType 1 Test: \n H0: all restrictions active(=)", 
            "\n     vs. \n H1: at least one restriction strictly true (>)", 
            "\n")
        print(aus.test, quote = FALSE)
        if (!brief) {
            cat("\nRestrictions on ", namen[object$restr.index[colSums(!object$ui == 
                0) > 0]], fill=TRUE)
            print(aus.rest, quote = FALSE, scientific = FALSE)
        }
        cat("\nRestricted estimate under H0:\n")
        print.default(format(object$b.eqrestr, digits = digits), 
            print.gap = 2, quote = FALSE)
        cat("\nRestricted estimate under union of H0 and H1 :\n")
        print.default(format(object$b.restr, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    if (object$TP == 11) {
        cat("\nType 11 Test: ", 
            "\n H0: all original restrictions active plus additional equality restrictions", 
            "\n     vs. \n H1: original restrictions hold", "\n")
        print(aus.test, quote = FALSE)
        if (!brief) {
            cat("\nRestrictions on", namen[object$restr.index[colSums(!object$ui == 
                0) > 0]], "\n")
            print(aus.rest, quote = FALSE, scientific = FALSE)
        }
        cat("\nRestricted estimate under union of H0 and H1 :\n")
        print.default(format(object$b.restr, digits = digits), 
            print.gap = 2, quote = FALSE)
        hilf <- object$ui.extra
        hilf[abs(hilf) < tol] <- 0
        colnames(hilf) <- namen
        rownames(hilf) <- paste("E", 1:nrow(hilf), ":", sep = "")
        hilf <- format(hilf, digits = digits)
        aus <- cbind(hilf, rep("%*%colnames == 0", nrow(object$ui.extra)))
        if (!brief) {
            cat("\nAdditional restrictions for H0:\n")
            print(aus, quote = FALSE)
        }
        cat("\nRestricted estimate under H0:\n")
        print.default(format(object$b.extra.restr, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    if (object$TP == 2) {
        cat("\nType 2 Test: \n H0: all restrictions true(>=)", 
            "\n     vs. \n H1: at least one restriction violated (<)", 
            "\n")
        print(aus.test, quote = FALSE)
        if (!brief) {
            cat("\nRestrictions on", namen[object$restr.index[colSums(!object$ui == 
                0) > 0]], "\n")
            print(aus.rest, quote = FALSE, scientific = FALSE)
        }
        cat("\nRestricted estimate under H0:\n")
        print.default(format(object$b.restr, digits = digits), 
            print.gap = 2, quote = FALSE)
        cat("\nUnrestricted estimate:\n")
        print.default(format(object$b.unrestr, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    if (object$TP == 21) {
        cat("\nType 21 Test: \n H0: all restrictions true(>= or =)", 
            "\n     vs. \n H1: at least one restriction violated (<), some =-restrictions maintained", 
            "\n")
        print(aus.test, quote = FALSE)
        if (!brief) {
            cat("\nRestrictions on", namen[object$restr.index[colSums(!object$ui == 
                0) > 0]], "\n")
            print(aus.rest, quote = FALSE, scientific = FALSE)
        }
        cat("\nRestricted estimate under H0:\n")
        print.default(format(object$b.restr, digits = digits), 
            print.gap = 2, quote = FALSE)
        cat("\nRestricted estimate under H1:\n")
        print.default(format(object$b.alt, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    if (object$TP == 3) {
        cat("\nType 3 Test: \n H0: at least one restriction not strictly true (<=)", 
            "\n         vs. \n H1: all restrictions strictly true (>)", 
            "\n")
        print(aus.test, quote = FALSE)
        if (!brief) {
            cat("\nRestrictions on", namen[object$restr.index[colSums(!object$ui == 
                0) > 0]], "\n")
            print(aus.rest, quote = FALSE, scientific = FALSE)
        }
        cat("\nUnrestricted estimate:\n")
        print.default(format(object$b.unrestr, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
}
