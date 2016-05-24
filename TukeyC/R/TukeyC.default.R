##
## S3 method to design matrix and response variable or data.frame objects
##

TukeyC.default <- function(x,
                           y=NULL,
                           model,
                           which,
                           error,
                           sig.level=.05,
                           round=2,
                           dispersion=c('mm', 's', 'se'), ...)   
{
  if (is.data.frame(y)) 
    y <- as.matrix(y[, 1])  # manova is not contemplated
  else
    stopifnot(is.atomic(y))

  if (is.matrix(x) || is.atomic(x))
    x <- as.data.frame(x)

  if(!is.null(y))
    dat <- as.data.frame(cbind(x,
                               y))
  else
    dat <- x

  av <- eval(substitute(aov(fo,
                            dat),
                        list(fo=formula(model))))

  if(class(av)[1] == 'aov')
    res <- TukeyC.aov(x=av,
                      which=which,
                      sig.level=sig.level,
                      round,
                      dispersion=dispersion)
  else
    res <- TukeyC.aovlist(x=av,
                          which=which, 
                          error=error,
                          sig.level=sig.level,
                          round,
                          dispersion=dispersion)

  invisible(res)
}
