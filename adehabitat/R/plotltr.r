"plotltr" <-
function(x, which="dist",...)
  {
    if (!inherits(x,"ltraj"))
      stop("x should be of class ltraj")
    opar <- par(mfrow=n2mfrow(length(x)))
    on.exit(par(opar))
    toto <- lapply(x, function(i) {
      ex<- parse(text=which)
      coin <- eval(ex, envir=i)
      plot(i$date, coin, main=attr(i,"burst"), xlab="Time",
           ylab=which, pch=16, cex=0.7,...)
      lines(i$date, coin)
    })
    invisible()
  }

