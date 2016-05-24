head.parametersMult <-
function(x, ...)
{
    for(iii in 1:x$dta.types)
        {
            x$dta[[iii]]$cg  <- head(x$dta[[iii]]$cg)
            x$dta[[iii]]$geo <- head(x$dta[[iii]]$geo)
        }
    print(x, ...)
}
