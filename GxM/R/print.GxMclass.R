print.GxMclass = function(GxMmle) {
    cat("Log-likelihood:\t", round(slot(GxMmle,"loglikelihood"),digits=2));
    cat("\n\n");

    cat("BIC:\t", round(slot(GxMmle,"BIC"),digits=2));
    cat("\n\n");

    cat("Coefficients:\n")
    print.default(round(slot(GxMmle,"par"), digits=4), print.gap=2, quote = FALSE)
    cat("\n")
    invisible(GxMmle);
}
