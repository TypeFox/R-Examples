tail.parametersSing <-
function(x, ...)
{
    x$dta$cg  <- tail(x$dta$cg)
    x$dta$geo <- tail(x$dta$geo)
    print(x, ...)
}
