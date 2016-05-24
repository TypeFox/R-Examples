geodistance <- function(longvar,latvar,lotarget,latarget,dcoor=FALSE) {

  latvar <- 2*pi*latvar/360
  longvar <- 2*pi*longvar/360
  lotarget <- 2*pi*lotarget/360
  latarget <- 2*pi*latarget/360
  dnorth <- NULL
  deast <- NULL
  dist <- pmin(sin(latvar)*sin(latarget) + cos(latvar)*cos(latarget)*cos(lotarget-longvar),  1)
  dist <- acos(dist)*3958

  if (dcoor==TRUE) {
    dnorth <- pmin(sin(latvar)*sin(latarget) + cos(latvar)*cos(latarget)*cos(0),  1)
    dnorth <- acos(dnorth)*3958
    dnorth <- ifelse(latvar<latarget, -dnorth, dnorth)
  }
  if (dcoor==TRUE) {
    deast <- pmin(sin(latvar)^2 + (cos(latvar)^2)*cos(lotarget-longvar),  1)
    deast <- acos(deast)*3958
    deast <- ifelse(longvar<lotarget, -deast, deast) 
  }
  out <- list(dist,dnorth,deast)
  names(out) <- c("dist","dnorth","deast")
  return(out)
}

