##' Conversion to long format data.frame for plotting
##'
##' Conversion to long format data.frame for plotting
##' @title dielectric2plot
##' @param m data.frame with wavelength and complex epsilon
##' @return long format data.frame
##' @author baptiste Auguie
##' @export
dielectric2plot <- function(m){

  dwide <- with(m, data.frame(wavelength,
                              real=Re(epsilon),
                              imag=Im(epsilon)))
  ## require(plyr)
  ## melt(dwide, id="wavelength")

  m <- 
    reshape(dwide, varying = c("real", "imag"),
            v.names="value", timevar="variable", 
            direction="long")
  m$variable <- factor(m$variable, labels=c("real", "imag"))
  m

}
