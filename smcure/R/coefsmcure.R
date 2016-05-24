##' @export
coefsmcure <-
function(x, ...)
{
	coef <- c(x$b,x$beta)
  names(coef) <- c(x$bnm,x$betanm)
  coef
}

