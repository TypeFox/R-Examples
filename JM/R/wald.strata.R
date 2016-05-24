wald.strata <-
function (fit) {
    if (!inherits(fit, "jointModel") || (nlv <- length(levels(fit$y$strata))) == 1) {
        stop("used only for stratified joint models.\n")
    }
    knots <- do.call(cbind, fit$control$knots)
    if (!all(apply(knots, 1, function (x) all(x == x[1])))) {
        warning("it appears that the knots are not the same among strata.\n")
    }
    coefs <- fit$coefficients
    spline.coefs <- coefs$gammas.bs
    ii <- grep("T.bs", colnames(fit$Hessian), fixed = TRUE)
    Var.spline.coefs <- vcov(fit)[ii, ii]
    n <- length(spline.coefs)
    p <- n / nlv
    L <- matrix(0, p * (nlv - 1), n)
    ind <- matrix(seq_len(n), ncol = nlv)
    pairs <- seq_len(nlv)
    pairs <- c(rbind(pairs[-nlv], pairs[-1]))    
    ii <- cbind(rep(seq_len(p), each = 2*(nlv - 1)), rep(pairs, p))
    L[cbind(rep(seq_len(nrow(L)), each = 2), ind[ii])] <- c(1, -1)
    v <- c(L %*% spline.coefs)
    stat <- c(crossprod(v, solve(L %*% tcrossprod(Var.spline.coefs, L), v)))
    df <- nrow(L)
    pvalue <- pchisq(stat, df = df, lower.tail = FALSE)
    mat <- rbind(c("X^2" = stat, df = df, "Pr(> X^2)" = pvalue))
    rownames(mat) <- rep("", nrow(mat))
    structure(list(alternative = "spline coefficients for the baseline risk\n\tfunction are not equal among strata", 
        Result = mat), class  = "wald.strata")
}
