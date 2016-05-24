tail.parametersCD <-
function(x, ...)
{
    x$dta.pc$cg <- tail(x$dta.pc$cg)
    x$dta.ph$cg <- tail(x$dta.ph$cg)
    x$dta.sc$cg <- tail(x$dta.sc$cg)
    x$dta.sh$cg <- tail(x$dta.sh$cg)

    x$dta.pc$geo <- tail(x$dta.pc$geo)
    x$dta.ph$geo <- tail(x$dta.ph$geo)
    x$dta.sc$geo <- tail(x$dta.sc$geo)
    x$dta.sh$geo <- tail(x$dta.sh$geo)

    print(x, ...)
}
