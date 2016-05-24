## -----------------------------------------------------------------------------
## Convert practical salinity to absolute salinity and vice versa
## -----------------------------------------------------------------------------

## internal functions
sw_delta_SA <- function (p0 = 0, lon0 = 0, lat0 = 0) {
  res <- NULL
  sw_sfac <- marelac::sw_sfac # import package data set
  for (i in 1:length(lon0))
    for (j in 1:length(lat0)){
      RES <-.Fortran("gsw_delta_sa", as.double(p0),
           as.double(lon0[i]), as.double(lat0[j]),
           sw_sfac$longs, sw_sfac$lats,
           sw_sfac$p, sw_sfac$ndepth, sw_sfac$del_sa,
           delta=as.double(0.))
       res <- c(res, RES$delta)
    }
    res
}


## Checks whether a lon, lat coordinate is in the Baltic
is_Baltic <- function(lon, lat) {
  leftx <- c(12.6, 7, 26)
  lefty <- c(50, 59, 69)

  rightx <- c(45, 26)
  righty <- c(50, 69)

  Baltic <- rep(FALSE, length(lat) * length(lon))
  ij <- 1
   for (i in 1:length(lon))
     for (j in 1:length(lat))  {
       if (lon[i] > min(leftx) & lon[i] < max(rightx) &
           lat[j]  > min(lefty) & lat[j]  < max(righty)) {
          xxl <- approx(lefty, leftx, lat[j])$y
          xxr <- approx(righty, rightx, lat[j])$y
          if (lon[i] > xxl & lon[i] < xxr) Baltic[ij] <- TRUE
       }
       ij <- ij+1
    }
  Baltic
}

## Salinity anomaly as a function of latitude, longitude, DSi conc and
## the ocean
DeltaSal <- function(p, lat, lon, DSi,
  Ocean = c("Global","Atlantic", "Pacific", "Indian", "Southern")) {

  dSal <- 0

  if (! is.null(DSi)) {  # a function of inputted DSi concentration
    Ocean <- match.arg(Ocean)
    if (Ocean == "Southern")
      dSal <-  7.4884e-5 *DSi
    else if (is.null(lat) | Ocean == "Global") {
      dSal <- 9.824e-5  *DSi
      if (Ocean != "Global")
        warning("latitude, lat, not given - using Global estimate instead")
      }
    else if (lat < -30)
      dSal <-  7.4884e-5 *DSi
    else if (Ocean == "Atlantic")
      dSal <- 7.4884e-5  *(1+1.0028*(lat/30+1))*DSi
    else if (Ocean == "Indian")
      dSal <- 7.4884e-5  *(1+0.3861*(lat/30+1))*DSi
    else if (Ocean == "Pacific")
      dSal <- 7.4884e-5  *(1+0.3622*(lat/30+1))*DSi

  } else if (! (is.null(lat)) | ! (is.null(lon))) {
    if (is.null(lat))
      stop (" latitude should be given if longitude is")
    if (is.null(lon))
      stop (" longitude should be given if latitude is")
    if (any(lat > 90 | lat < -90))
      stop (" latitude, 'lat' should be within -90,+90")

    if (any(lon > 360 | lon < 0))
      stop (" longitude, 'lon' should be within 0,360")

    dSal <- sw_delta_SA(p * 10, lon, lat)

  }

  dSal
}

## Practical -> absolute salinity
convert_PStoAS <- function(S=35, p=max(0,P-1.013253), P=1.013253,
  lat=NULL, lon=NULL, DSi = NULL,
  Ocean = c("Global", "Atlantic", "Pacific", "Indian", "Southern")) {

  SA   <- 35.16504 /35.0 *S +DeltaSal(p, lat, lon, DSi, Ocean)

  if (is.null(DSi) & ! (is.null(lat)) & ! (is.null(lon))) {
    Baltic <- is_Baltic(lon, lat)
    if (length(Baltic)!=length(SA)) Baltic <- rep(Baltic, len=length(SA))
    SA[Baltic] <- SA[Baltic] + 0.124*(1-S/35)
  }
  SA
}

## Absolute -> practical salinity
convert_AStoPS <- function(S = 35, p = max(0, P-1.013253), P = 1.013253,
  lat = NULL, lon = NULL, DSi = NULL,
  Ocean = c("Global", "Atlantic", "Pacific", "Indian", "Southern")) {

  PS <- 35.0/35.16504 *(S - DeltaSal(p, lat, lon, DSi, Ocean))

  if (is.null(DSi) & ! (is.null(lat)) & ! (is.null(lon))) {
    Baltic <- is_Baltic(lon, lat)
    if (length(Baltic) != length(PS)) Baltic <- rep(Baltic, len = length(PS))
    PS[Baltic] <- PS[Baltic] - 35.0/35.16504 * S + (35.0/35.04104) * (S - 0.124)
  }
  PS
}

