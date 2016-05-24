`is.fts` <- function (x) 
{
    inherits(x, "fts") & length(x$x) > 0
}
