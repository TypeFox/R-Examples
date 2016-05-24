print.rcor.test <-
function (x, digits = max(3, getOption("digits") - 4), ...) {
    mat <- x$cor.mat
    mat[lower.tri(mat)] <- x$p.values[, 3]
    mat[upper.tri(mat)] <- sprintf("%6.3f", as.numeric(mat[upper.tri(mat)]))
    mat[lower.tri(mat)] <- sprintf("%6.3f", as.numeric(mat[lower.tri(mat)]))
    ind <- mat[lower.tri(mat)] == paste(" 0.", paste(rep(0, digits), collapse = ""), sep = "")
    mat[lower.tri(mat)][ind] <- "<0.001"
    ind <- mat[lower.tri(mat)] == paste(" 1.", paste(rep(0, digits), collapse = ""), sep = "")
    mat[lower.tri(mat)][ind] <- ">0.999"
    diag(mat) <- " *****"
    cat("\n")
    print(noquote(mat))
    cat("\nupper diagonal part contains correlation coefficient estimates",
        "\nlower diagonal part contains corresponding p-values\n\n")
    invisible(x)
}
