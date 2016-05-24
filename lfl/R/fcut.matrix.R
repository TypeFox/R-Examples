fcut.matrix <- function(x, ...) {
    result <- fcut(as.data.frame(x), ...)
    return(result)
}
