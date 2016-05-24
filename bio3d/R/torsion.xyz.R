"torsion.xyz" <-
function(xyz, atm.inc=4) {

  if(!is.vector(xyz) || !is.numeric(xyz))
    stop("input 'xyz' should be a numeric vector") 
  natm  <- length(xyz)/3
  if(natm < 4)
    stop("Need at least four atoms to define a dihedral")
  if(natm %% 1 != 0)
    stop("There should be three 'xyz' elements per atom")

  m.xyz <- matrix(xyz, nrow=3)
  atm.inds <- c(1:4); out<-NULL
  while(atm.inds[4] <= natm) {
    if( any(is.na( m.xyz[,atm.inds] )) ) {
      torp <- NA
    } else {
      d1 <- m.xyz[,atm.inds[2]] - m.xyz[,atm.inds[1]]
      d2 <- m.xyz[,atm.inds[3]] - m.xyz[,atm.inds[2]]
      d3 <- m.xyz[,atm.inds[4]] - m.xyz[,atm.inds[3]]

      u1 <- (d1[c(2,3,1)] * d2[c(3,1,2)]) - (d2[c(2,3,1)] * d1[c(3,1,2)])
      u2 <- (d2[c(2,3,1)] * d3[c(3,1,2)]) - (d3[c(2,3,1)] * d2[c(3,1,2)])

      ctor <- sum(u1*u2)/sqrt( sum(u1*u1) * sum(u2*u2) )
      ctor[ctor > 1] <- 1; ctor[ctor < -1] <- -1
      torp <- matrix(acos(ctor)*(180/pi),ncol=1)

      if( sum(u1 * ((u2[c(2,3,1)] * d2[c(3,1,2)]) -
                    (u2[c(3,1,2)] * d2[c(2,3,1)]))) < 0)
        torp <- -torp
    }
    out <- c(out, torp)
    atm.inds <- atm.inds + atm.inc
  }
  if(atm.inc == 1 & natm > 4) out <- c(NA, out, NA, NA)
  return(out)
}

