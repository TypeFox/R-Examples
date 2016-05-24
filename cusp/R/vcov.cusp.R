`vcov.cusp` <-
function (object, ...) 
{
    so <- summary.cusp(object, correlation = FALSE, logist = FALSE, 
        ...)
    so$dispersion * so$cov.unscaled
}

