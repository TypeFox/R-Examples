"angle.xyz" <-
function(xyz, atm.inc=3) {
  if(!is.vector(xyz) || !is.numeric(xyz))
    stop("input 'xyz' should be a numeric vector")
  natm  <- length(xyz)/3
  if(natm < 3)
    stop("Need at least three atoms to define an angle")
  if(natm %% 1 != 0)
    stop("There should be three 'xyz' elements per atom")

  m.xyz <- matrix(xyz, nrow=3)
  atm.inds <- c(1:3); out <- NULL
  while(atm.inds[3] <= natm) {
    if( any(is.na( m.xyz[,atm.inds] )) ) {
      ang <- NA
    } else {
      d1 <- m.xyz[,atm.inds[1]] - m.xyz[,atm.inds[2]]
      d2 <- m.xyz[,atm.inds[3]] - m.xyz[,atm.inds[2]]
      ang <- sum(d1*d2) / (sqrt(sum(d1^2)) * sqrt(sum(d2^2)) )
      ang[ang > 1] <- 1; ang[ang < -1] <- -1
      ang <- acos(ang) * (180/pi)
    }
    out <- c(out, ang)
    atm.inds <- atm.inds + atm.inc
  }
  return(out)
}

