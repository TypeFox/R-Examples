print.magpie <- function(x, drop=TRUE, ...) {
    print(as.array(x)[,,,drop=drop], ...)
}