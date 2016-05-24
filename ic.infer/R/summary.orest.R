summary.orest <- function (object, display.unrestr = FALSE, brief = FALSE, digits = max(3, 
    getOption("digits") - 3), scientific = FALSE, ...) 
{
    cat("Order-restricted estimated mean with restrictions of coefficients of", 
        "\n")
    cat(names(object$b.restr)[object$restr.index[colSums(!object$ui == 0) > 
        0]], "\n\n")
    namen <- names(object$b.restr)
    if (is.null(namen)) 
        namen <- paste("m", 1:length(object$b.restr), sep = "")
    hilf <- object$ui
    colnames(hilf) <- namen[object$restr.index]
    aus <- cbind(hilf, rep("%*%colnames", nrow(object$ui)), c(rep(" == ", 
        object$meq), rep(" >= ", nrow(object$ui) - object$meq)), object$ci)
    rownames(aus) <- paste(1:nrow(aus), ":", sep = "")
    colnames(aus)[(ncol(hilf) + 1):ncol(aus)] <- rep(" ", ncol(aus) - 
        ncol(hilf))
    aus <- cbind(rep(" ", nrow(aus)), aus)
    aus[object$iact, 1] <- "A"
    if (!brief) {
        if (object$meq == 0) 
            cat("\n", "Inequality restrictions:\n")
        else cat("\n", "Restrictions:\n")
        print(aus, quote = FALSE, scientific = FALSE)
        cat("\n", "Note: Restrictions marked with A are active.", 
            "\n")
    }
    g <- length(object$b.restr)
    cat("\nRestricted estimate:", "\n")
    orc <- round(object$b.restr, 9)
    mark <- rep(" ", length(object$b.restr))
    mark[object$restr.index[colSums(!object$ui == 0) > 0]] <- "R"
    names(orc) <- paste(mark, names(orc), sep = " ")
    print(format(orc), quote = FALSE, digits = 4)
    cat("\n", "Note: Estimates marked with R are involved in restrictions.", 
        "\n")
    if (display.unrestr) {
        cat("\n\nUnrestricted estimate: ", "\n")
        print(object$b.unrestr)
    }
}
