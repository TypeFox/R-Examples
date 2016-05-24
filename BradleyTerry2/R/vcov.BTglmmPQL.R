vcov.BTglmmPQL <- function (object, ...)
{
    so <- summary(object, corr = FALSE, ...)
    so$dispersion * so$cov.unscaled
}
