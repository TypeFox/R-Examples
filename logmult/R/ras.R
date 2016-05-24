ras <- function(tab, row, col, tolerance=.Machine$double.eps) {
    if(!isTRUE(all.equal(sum(row), sum(col))))
        stop("sum(u) must be equal to sum(v)")

    if(any(tab == 0))
        warning("convergence is not guaranteed when some cells are equal to 0")

    if(any(tab < 0))
        stop("elements of tab must all be >= 0")

    if(any(row <= 0) || any(col <= 0))
        stop("elements of row and col must all be > 0")

    # Destroyed by operations
    attr <- attributes(tab)

    while(!(isTRUE(all.equal(rowSums(tab), row, tolerance)) &&
            isTRUE(all.equal(colSums(tab), col, tolerance)))) {
        r <- row / rowSums(tab)
        tab <- diag(r) %*% tab
        c <- col / colSums(tab)
        tab <- tab %*% diag(c)
    }

    attributes(tab) <- attr
    tab
}
