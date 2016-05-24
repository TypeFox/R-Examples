##
## S3 method to 'aovlist' object
##

TukeyC.aovlist <- function(x,
                           which,
                           error,
                           sig.level=.05,
                           round=2,
                           dispersion=c('mm', 's', 'se'), ...)
{
  mt <- model.tables(x,
                     "means")              # summary tables for model fits
  if(is.null(mt$n))
    stop("No factors in the fitted model!")

  tabs <- mt$tables[-1][which]             # specified group means

  r    <- mt$n[names(tabs)][[which]]       # groups and its number of replicates

  bal  <- ifelse(length(r) == 1,
                 TRUE,
                 FALSE)                    # is (or not) balanced

  MSE  <- sum(resid(x[[error]])^2) / x[[error]][[8]]

  nms  <- names(tabs[[which]])

  ord  <- order(as.vector(tabs[[which]]),
                decreasing=TRUE)

  m.inf <- m.inf.1b(x,
                    which,
                    dispersion)
  
  rownames(m.inf) <- nms

  m.inf <- m.inf[order(m.inf[,1],
                       decreasing=TRUE),]

  dfr <- x[[error]][[8]]  # residual degrees of freedom

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

