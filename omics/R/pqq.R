pqq <- function(ps) {
    ps <- na.omit(ps)
    n <- length(ps)
    list(
        x=-log10(1:n / n),
        y=-log10(sort(ps))
    )
}

pqq.ci <- function(n, level=0.95) {
    alpha <- (1 - level) / 2
    result <- -log10(cbind(
        qbeta(alpha, 1:n, n - 1:n + 1),
        qbeta(1 - alpha, 1:n, n - 1:n + 1)
    ))
    colnames(result) <- format.percentage(c(alpha, 1 - alpha), 3)
    result
}
