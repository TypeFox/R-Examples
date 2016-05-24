UTMzone<-function(lat, lon=NA)
  {
### return the utm zone and origin for projection

###   this is how the matla code works
    ##  there are two cases: if lat is a character convert UTM to lat lon
    ##  if lat is numeric: convert lat-lon to UTMzone

###  setting things up:
    utmx = seq(from = (-180), to = 180, by = 6)
    utmy = seq(from = -56, to = 72, by = 8)
    LONS = RPMG::fmod(utmx, 360)
    htags = LETTERS[6:23]
    htags = htags[-c(4, 10)]


####  case 1
    if(is.na(lon))
      {
        Kzone = lat
        K1 = unlist( strsplit(Kzone, split="") )
        klen = length(K1)
        h1 = K1[klen]
        h2 = as.numeric(paste(K1[1:(klen-1)], collapse=""))
        ILON = h2
        ILAT = which(h1==htags)
        latzone =  htags[ILAT]
        hh = paste( ILON, latzone, sep="")

        LONs   = c(utmx[ILON],utmx[ILON]+6)
        LONcen = mean(LONs)
        LATs   = c(utmy[ILAT], utmy[ILAT]+8)
        LATcen = mean( LATs)
        return(list(zone=hh, LON=LONs, LAT=LATs, CEN=c(LATcen, LONcen) ))
      }
    else
      {
        ILON   = findInterval( lon,  utmx  )

        ILAT = findInterval(lat,  utmy  )

        latzone =  htags[ILAT]

        hh = paste( ILON, latzone, sep="")

        
        LONs   = c(utmx[ILON],utmx[ILON]+6)
        LONcen = mean(LONs)
        LATs   = c(utmy[ILAT], utmy[ILAT]+8)
        LATcen = mean( LATs)

        return(list(zone=hh, LON=LONs, LAT=LATs, CEN=c(LATcen, LONcen) ))
      }


  }

