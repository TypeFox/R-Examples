head.parametersSing <-
function(x, ...)
{
    x$dta$cg  <- head(x$dta$cg)
    x$dta$geo <- head(x$dta$geo)
    print(x, ...)
}
