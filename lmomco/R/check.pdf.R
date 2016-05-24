"check.pdf" <-
function(pdf, para, lowerF=0.001, upperF=0.999, eps=0.02, verbose=FALSE,
         plot=FALSE, plotlowerF=0.001, plotupperF=0.999, ...) {

  if(lowerF > upperF) {
    tmp <- lowerF
    lowerF <- upperF
    upperF <- tmp
    warning("End points appear backwards, swapping for you")
  }

  if(lowerF == upperF) {
    warning("End points are same, returning NULL")
    return(NULL)
  }

  min <- par2qua(lowerF,para)
  max <- par2qua(upperF,para)

  if(max < min) {
    warning("quantile function max < min, returning NULL")
    return(NULL)
  }

  # Perform the integration using defaults
  pdfgral <- integrate(pdf,min,max,para=para,...)

  if(verbose) {
    cat(c("STATUS: integrated to ", pdfgral$value,
          " with absolute error <", pdfgral$abs.error,"\n"))
  }

  unity <- FALSE

  if(abs(pdfgral$value-1) < eps) {
     print("pdf function appears to integrate to unity")
     unity <- TRUE
  }
  else {
     warning("pdf function does not integrate to unity. ",
             "You might considered plotting the function and looking for far ",
             "tail performance problems of the function or you otherwise have ",
             "poor end points for the numerical integration")
  }

  if(plot) {
     subs <- (plotupperF - plotlowerF)/1000
     F <- seq(plotlowerF,plotupperF,by=subs)
     F <- F[F >= 0]
     F <- F[F <= 1]
     x <- par2qua(F,para)
     y <- pdf(x,para); y[is.na(y)] <- 0 # and reset NAs to zero
     plot(x,y,type='l',ylab="probability density")

     min1 <- min(y)

     subs <- (upperF - lowerF)/1000
     F <- seq(lowerF,upperF,by=subs)
     F <- F[F >= 0]
     F <- F[F <= 1]
     x <- par2qua(F,para)
     y <- pdf(x,para); y[is.na(y)] <- 0 # and reset NAs to zero

     min2 <- min(y)
     min <- min(min1,min2)

     polygon(c(x[1],x,max(x)),c(min,y,min), col=rgb(.5,.5,.5))
  }

  return(list(isunity=unity, F=pdfgral$value))
}
