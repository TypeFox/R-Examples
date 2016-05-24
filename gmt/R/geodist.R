geodist <- function(Nfrom, Efrom, Nto, Eto, units="km")
{
  units <- match.arg(units, c("km","nm"))

  rad <- 180 / pi

  N1 <- Nfrom / rad
  E1 <- Efrom / rad
  N2 <- Nto   / rad
  E2 <- Eto   / rad

  duplicates <- N1==N2 & E1==E2
  N1[duplicates] <- 0            # When origin and destination are the same,
  E1[duplicates] <- 0            #   set them both to 0, 0
  N2[duplicates] <- 0            # Without this, geodist(48.535, 124, 48.535, 124) returns NaN,
  E2[duplicates] <- 0            #   but geodist(0, 0, 0, 0) seems to return 0 on all machines

  radians <- acos(sin(N1)*sin(N2)+cos(N1)*cos(N2)*cos(E1-E2))

  distance <- if(units=="km") 60*rad*radians*1.852 else 60*rad*radians

  return(distance)
}
