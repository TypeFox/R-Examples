##
## S3 method to design matrix and response variable or data.frame objects
##

SK.nest.default <- function(x,
                            y=NULL,
                            model,
                            which,
                            id.trim=3,
                            error,
                            fl1,
                            fl2=0,
                            sig.level=.05,
                            dispersion=c('mm', 's', 'se'), ...)
{ 
  if (is.data.frame(y))
    y <- as.matrix(y[, 1])  # manova is not contemplated
  else
    stopifnot(is.atomic(y))

  if (is.matrix(x) || is.atomic(x))
    x <- as.data.frame(x)

  if(!is.null(y))
    dat <- as.data.frame(cbind(x, y))
  else
    dat <- x

  av <- eval(substitute(aov(fo,
                            x),
                        list(fo=formula(model))))

  if(class(av)[1] == 'aov')
    res <- SK.nest.aov(x=av,
                       which=which,
                       id.trim=id.trim,
                       fl1=fl1,
                       fl2=fl2,
                       sig.level=sig.level,
                       dispersion=dispersion)
  else
    res <- SK.nest.aovlist(x=av,
                           which=which,
                           id.trim=id.trim,
                           error=error,
                           fl1=fl1,
                           fl2=fl2,
                           sig.level=sig.level,
                           dispersion=dispersion)

  invisible(res)
}
