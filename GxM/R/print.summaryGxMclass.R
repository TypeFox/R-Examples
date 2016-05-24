print.summaryGxMclass <-
function(x, ..., eps=1e-10) {

    sumGxMmle = x;
    cat("Coefficients:\n")
    est = slot(sumGxMmle,"par");
    estnames = names(est);
    stderr = slot(sumGxMmle,"stderr");
    zvalues = ifelse(stderr > eps, est/stderr, NA);
    cmat = cbind(est, stderr, zvalues, 2*pnorm(-abs(zvalues)));
    colnames(cmat) = c("Estimate", "Std Err", "z value", "Pr(>|z|)")
    rownames(cmat) = names(est);
    printCoefmat(cmat);
    cat("\n");
    cat("Log-likelihood:\t", round(slot(sumGxMmle,"loglikelihood"),digits=2));
    cat("\n");
    cat("BIC:\t", round(slot(sumGxMmle,"BIC"),digits=2));
    cat("\n");
    invisible(sumGxMmle);
}
