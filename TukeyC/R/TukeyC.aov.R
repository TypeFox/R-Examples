##
## S3 method to 'aov' object
##

TukeyC.aov <- function(x,
                       which=NULL,
                       sig.level=.05,
                       round=2,
                       dispersion=c('mm', 's', 'se'), ...)   
{
  if(is.null(which))
    which <- names(x$model)[2]

  mt <- model.tables(x,
                     "means")               # summary tables for model fits

  if(is.null(mt$n))
    stop("No factors in the fitted model!")

  tabs <- mt$tables[-1][which]              # specified group means                             

  r    <- mt$n[names(tabs)][[which]]        # groups and its number of replicates

  bal  <- ifelse(length(r) == 1,
                 TRUE,
                 FALSE)                     # is (or not) balanced

  MSE  <- sum(resid(x)^2) / x$df.residual

  nms  <- names(tabs[[which]])

  ord  <- order(as.vector(tabs[[which]]),
                decreasing=TRUE)

  m.inf <- m.inf.1a(x,
                    which,
                    dispersion)

  rownames(m.inf) <- nms  

  m.inf <- m.inf[order(m.inf[,1],
                       decreasing=TRUE),]

  dfr <- x$df.residual  # residual degrees of freedom

  out <- make.TukeyC.test(r=r,
                          MSE=MSE,
                          m.inf=m.inf,
                          ord=ord,
                          sig.level=sig.level,
                          dfr=dfr,
                          bal=bal,
                          mt=mt,
                          round)

  class(out) <- c('TukeyC',
                  'list')

  return(out) 
}
