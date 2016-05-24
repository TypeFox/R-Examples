##
## S3 method to design matrix and response variable or data.frame objects
##

TukeyC.nest.default <- function(x,
                                y=NULL,
                                model,
                                which,
                                error,
                                fl1,
                                fl2=0,
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
                            x),
                        list(fo=formula(model))))

  if(class(av)[1] == 'aov')
    res <- TukeyC.nest.aov(x=av,
                           which=which,
                           fl1=fl1,
                           fl2=fl2,
                           sig.level=sig.level,
                           round)
  else
    res <- TukeyC.nest.aovlist(x=av,
                               which=which,
                               error=error,
                               fl1=fl1,
                               fl2=fl2,
                               sig.level=sig.level,
                               round)
  invisible(res)
}
