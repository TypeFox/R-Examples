PWV <- function(sonde=NULL, press, td, temp, minp = 400.)
{
# calculate precipitable water (in mm) up to minp (minimum pressure) using
# dew-point temperature (td) and temperature (temp) (both in deg C)
#
# Grab data from sonde if not null, otherwise use values passed in.
if( !is.null( sonde) ) {
	press <- sonde$press
	td <- sonde$dewpt
	temp <- sonde$temp
	}

  q <- td.to.q(temp, td, press)
  n <- length(press)
  pw <- 0.
  for(k in 2.:n) {
    if(press[k] < minp)
      break
    if(is.na(q[k]) || is.na(q[k - 1.]))
      next
    meanq <- (q[k] + q[k - 1.])/2.
    pw <- pw + (0.1 * meanq * (press[k - 1.] - press[k]))/9.81
  }
  return(pw)
}
