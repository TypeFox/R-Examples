rbeta <- function (n, shape1, shape2) {
    if (length(n)>1) n <- length(n)
    .Call("newrbeta", n, shape1, shape2)
}
