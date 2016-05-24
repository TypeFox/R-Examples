Hubble.LCDM <- function(z,Omega.M0=0.3,Omega.L0=0.7){ sqrt( Omega.M0*(1.+z)^3+Omega.L0+(1.-Omega.M0-Omega.L0)*(1.+z)^2 ) }
