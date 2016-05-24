tail.parametersMult <-
function(x, ...)
{
    for(iii in 1:x$dta.types)
        {
            x$dta[[iii]]$cg  <- tail(x$dta[[iii]]$cg)
            x$dta[[iii]]$geo <- tail(x$dta[[iii]]$geo)
        }
    print(x, ...)
}
