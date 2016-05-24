print.mvformula <-
function (x, ...) 
{
    cat("mvformula\n")
    attr(x, ".Environment") <- NULL
    print.default(unclass(x), ...)
}
