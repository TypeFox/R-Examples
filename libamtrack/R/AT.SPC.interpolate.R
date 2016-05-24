AT.SPC.interpolate <- function(spc.lower, spc.upper, energy.MeV.u){
  if(spc.lower$n.depth.steps != spc.upper$n.depth.steps){
    # TODO: rescale spc with more steps to steps from spc with
    # TODO: fewer steps using AT.SPC.spectrum.at.depth.g.cm2
    cat("Cannot interpolate, different number of depth steps")
    return(NULL)
  }
  if((spc.lower$projectile != spc.upper$projectile)|
     (spc.lower$target.material != spc.upper$target.material)){
    cat("Cannot interpolate, target materials or projectiles do not match.")
    return(NULL)
  }
  if((spc.lower$energy.MeV.u > energy.MeV.u)|
     (spc.upper$energy.MeV.u < energy.MeV.u)){
    cat("Cannot interpolate, requested energy outside of given spc energies.")
    return(NULL)
  }
  # Copy structure
  spc                    <- spc.lower
  spc$energy.MeV.u       <- energy.MeV.u
  
  # Get relative difference of requested energy
  frac                   <- (energy.MeV.u - spc.lower$energy.MeV.u) / (spc.upper$energy.MeV.u - spc.lower$energy.MeV.u)
  
  # scale depth steps
  # TODO: should scales with E2 like range
  spc$spc$depth.g.cm2    <- (1-frac) * spc.lower$spc$depth.g.cm2 + frac * spc.upper$spc$depth.g.cm2 

  # interpolate fluences
  # TODO: check if linear interpolation really applies
  spc$spc$fluence.cm2    <- (1-frac) * spc.lower$spc$fluence.cm2 + frac * spc.upper$spc$fluence.cm2

  return(spc)
}