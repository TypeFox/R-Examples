"[.snp"<-
function (x, ..., drop = FALSE)
{
    y <- NextMethod("[")
#    attr(y, "contrasts") <- attr(x, "contrasts")
#    attr(y, "levels") <- attr(x, "levels")
#    class(y) <- oldClass(x)
    attributes(y) <- attributes(x)
    if (drop)
        factor(y)
    else y
}
