# Computation of the ARL for modified Shewhart charts, AR(1) data
xshewhart.ar1.arl <- function(alpha, cS, delta=0, N1=50, N2=30) {
  if ( abs(alpha) >= 1 )    stop("alpha has to be between -1 and 1")
  if ( cS <= 0 )            stop("cS has to be positive")
  arl <- .C("xshewhart_ar1_arl",
              as.double(alpha), as.double(cS), as.double(delta),
              as.integer(N1), as.integer(N2), ans=double(length=1), PACKAGE="spc")$ans
  names(arl) <- NULL
  arl
}
