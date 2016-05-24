"itCall" <- function(symbol, ...)
{
    args <- list(...)
    CLASSES <- as.character(sapply(args, function(x) class(x)[1L]))
    COPY <- rep(FALSE, length(args))
    .Call(symbol, ..., COPY = COPY, CLASSES = CLASSES, PACKAGE = "ifultools")
}
