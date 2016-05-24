print.SPRT <-
function(x = SPRT,...) {
    cat("Wald's Sequential Probability Ratio Test (SPRT)\n\n")
    
    if (!is.null(x$decision)) {
        cat("Decision:", x$interpretation, "", sep = "\n")
        cat("Distribution:", x$distribution, "\n", sep = " ")
        cat("n: ", x$n, ", k: ", x$k, "\n", sep = "")
        cat("h0: ", x$h0, ", h1: ", x$h1, "\n\n", sep = "")
    }

    cat("Wald boundaries (log):\n")
    cat("> B boundary:      ", round(x$wald.B, 3), "\n", sep = " ")
    cat("> A boundary:       ", round(x$wald.A, 3), "\n", sep = " ")
    if (!is.null(x$llr)) {
        cat("> Likelihood ratio: ", round(x$llr, 3), "\n", sep = " ")   
    }
    if (!is.null(x$data.sum)) {
        cat("\nPreview k boundaries:\n")
        print.data.frame(head(x$data.sum, 5), digits = 3, row.names = FALSE) 
    }
}
