`summary.roc` <- function(object, ...) {
    summ <- object$statistics
    class(summ) <- "summary.roc"
    return(summ)
}
