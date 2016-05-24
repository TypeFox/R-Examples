digits <- function(x, n=NULL, simplify = FALSE) {
    if(length(x) > 1) {
        if(is.null(n) & simplify) {
            n <- floor(max(log10(x))) + 1
        }
        sapply(x,digits, simplify=simplify, n=n)
    } else {
        if(is.null(n)) {
            n <- floor(log10(x))
        } else {
            n <- n - 1
        }
        x %/% 10^seq(n,0) %% 10
    }
}
