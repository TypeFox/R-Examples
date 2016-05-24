summary.dlc <- function(object, ...){

    x <- object$x
    phi <- object$phi
    n <- length(object$xn)
    m <- length(x)

    mode <- x[phi == max(phi)]
    phi.mode <- phi[x == mode]
    f.mode <- exp(phi.mode)
    knots <- (1:m)[object$IsKnot == 1]

    p <- "\nEstimation of a log-concave density from i.i.d. data\n"

    pn <- paste("\nNumber of initial observations: n = ", n, "\n", sep = "")

    p0 <- paste("\nNumber of unique observations (or grid points): m = ", m, "\n", sep = "")

    pw <- paste("\nMinimal and maximal weight: ", min(object$w), ", ", max(object$w), "\n", sep = "")

    p1 <- paste("\nlog-likelihood: ", round(object$L, 2), "\n", sep = "")

    p2 <- paste("\nMaximum likelihood estimate:\nMode: x[", (1:m)[x == mode], "] = ", format(mode, digits = 2), "\n", "Value of log-density at mode: ", format(phi.mode, digits = 2),
        "\n", "Value of density at mode: ", format(f.mode, digits = 2), "\n", sep = "")

    p3 <- paste("\nNumber of knots of the MLE: ", sum(object$IsKnot), "\n", sep = "")

    p4 <- paste("\nKnots of the MLE:\nx[", paste(knots, collapse = ", "), "] =\n", paste(format(x[object$IsKnot == 1], digits = 2), collapse = ", "), "\n\n", sep = "")

    if (object$smoothed == FALSE){cat(p, pn, p0, pw, p1, p2, p3, p4)}

    if (object$smoothed == TRUE){
        mode.s <- object$xs[object$f.smoothed == max(object$f.smoothed)]
        f.mode.s <- object$f.smoothed[object$xs == mode.s]
        phi.mode.s <- log(f.mode.s)

        p5 <- paste("\nSmoothed maximum likelihood estimate:\nMode: ", format(mode.s, digits = 2), "\n",
            "Value of log-density at mode: ", format(phi.mode.s, digits = 2),
        "\n", "Value of density at mode: ", format(f.mode.s, digits = 2), "\n", sep = "")
        cat(p, pn, p0, p1, p2, p5, p3, p4)}       
}
