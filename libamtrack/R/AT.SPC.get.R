AT.SPC.get <- function(spc.list, energy.MeV.u){
  if(energy.MeV.u < min(spc.list$energy.MeV.u) | energy.MeV.u > max(spc.list$energy.MeV.u)){
    cat("Requested energy outside of grid of available data.")
    return(NULL)
  }
  
  # find two closest matches
  distance    <- abs(spc.list$energy.MeV.u - energy.MeV.u)
  closest.idx <- match(head(sort(distance), 2), distance)
  if(closest.idx[1] == closest.idx[2]){
    closest.idx[2] <- closest.idx[1] + 1
  }
  
  # read and interpolate SPCs
  spc.lower   <- AT.SPC.read( file.name        = spc.list$file.name[closest.idx[1]],,
                              endian           = spc.list$endian[closest.idx[1]])
  spc.upper  <- AT.SPC.read( file.name        = spc.list$file.name[closest.idx[2]],,
                              endian           = spc.list$endian[closest.idx[2]])

  spc         <- AT.SPC.interpolate(spc.lower    = spc.lower,
                                    spc.upper   = spc.upper,
                                    energy.MeV.u = energy.MeV.u)
  return(spc)
}
