"compModelMass" <-
function (theta, model)
{
  peakpar <- theta@peakpar
  amp <- theta@amplitudes
  if(model@extracomp) {
    if(length(amp)==0)
      eamp <- 1
    else {
      eamp <- amp[length(amp)]
      amp <- amp[- length(amp)]
    }
  }
  shift <- theta@shift
  if(model@peakfunct == "expmodgaus") {
    fn1 <- function(x,ind) x[[ind]]
    lpp <- length(peakpar)
    locations <- unlist(lapply(peakpar, fn1, ind=1))
    if(length(shift)!=0)
      locations <- locations + shift 
    widths <- unlist(lapply(peakpar, fn1, ind=2))
    rates <- unlist(lapply(peakpar, fn1, ind=3))
    massm <- rep(0, model@nt * lpp) 
    massm <- as.matrix(.C("calcCirf_multi", 
                          cmat = as.double(massm), 
                          as.double(rates), as.double(model@x), 
                          as.double(widths), 
                          as.double(locations), 
                          as.integer(lpp), 
                          as.integer(model@nt), PACKAGE="TIMP")$cmat)
    dim(massm) <- c(model@nt, lpp)
    if(length(amp) > 0)
      massm <- massm %*% diag(amp)
  }
  if(model@extracomp)
    massm <- cbind(massm, rep(eamp, model@nt))
  massm
}
