su <-
  function(days,lat,lon,extraT=NULL,A=NA,B=NA,C=NA,tmax,tmin,CC) {
    i <- dayOfYear(days)
    latt <- radians(lat)
    if (is.null(extraT)) extraT <- extrat(i=i,lat=latt)$ExtraTerrestrialSolarRadiationDaily
    
    if (is.na(A)) A <- extract(Sa_map, matrix(c(lon,lat),1,2))
    if (is.na(A)) stop("Lat/lon outside the coefficient map!")
    if (is.na(B)) B <- extract(Sb_map, matrix(c(lon,lat),1,2))
    if (is.na(B)) stop("Lat/lon outside the coefficient map!")
    if (is.na(C)) C <- extract(Sc_map, matrix(c(lon,lat),1,2))
    if (is.na(C)) stop("Lat/lon outside the coefficient map!")
    
    su <- extraT*(A*sqrt(tmax-tmin)+B*sqrt(1-CC/8)) + C
    su 
  }

