kalsmoComp <- 
  function (x, comp=NULL, plot.it = TRUE, xlab="time",ylab="",na.action = na.fail, ...) {
    cn <- match(c("sser"), names(x))
    if (any(is.na(cn))) 
      stop("x must be a kalsmo.car() object")
    if(length(comp)==1)
    compser <- Re(x$sser[,comp])
    else
    compser <- Re(apply(x$sser[,comp],1,sum))
    ntim <- length(x$ser)
    spg.out <- list(tim=x$tim[1:ntim],compser=compser, method = paste("CAR (", 
                               comp, ") component ", sep = ""))
    class(spg.out) <- "kalsmo.comp"
    if (plot.it) {
      plot(spg.out$tim,spg.out$compser,xlab=xlab,ylab=ylab,type="l",...)
      return(invisible(spg.out))
    }
    else return(spg.out)
  }
