head.parametersCD <-
function(x, ...)
{
    x$dta.pc$cg <- head(x$dta.pc$cg)
    x$dta.ph$cg <- head(x$dta.ph$cg)
    x$dta.sc$cg <- head(x$dta.sc$cg)
    x$dta.sh$cg <- head(x$dta.sh$cg)

    x$dta.pc$geo <- head(x$dta.pc$geo)
    x$dta.ph$geo <- head(x$dta.ph$geo)
    x$dta.sc$geo <- head(x$dta.sc$geo)
    x$dta.sh$geo <- head(x$dta.sh$geo)

    print(x, ...)
}
