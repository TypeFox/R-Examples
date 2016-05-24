"plot.ahaz"<-function(x, ...)
  {
    if(x$univar)
      stop("not supported for univariate regressions")
    out <- predict(x,type="cumhaz")
    plot(out, ...)
  }
