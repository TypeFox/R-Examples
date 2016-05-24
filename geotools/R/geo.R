##
## return approximated distance in kms
##

distKmRad<- function(lat0,lon0,lat1,lon1)
{
  deltalat = lat1 - lat0
  deltalon = lon1-lon0
  ## Earth radius (km)
  RT = 6400
  R = 6400 * cos(lat0)
  return( sqrt( (RT*deltalat)^2 + (R*deltalon)^2) )
}

distKm<- function(lat0,lon0,lat1,lon1)
{
  ## convert to radian
  lat0 = pi*lat0 / 180.0
  lat1 = pi*lat1 / 180.0
  lon0 = pi*lon0 / 180.0
  lon1 = pi*lon1 / 180.0
  return (distKmRad(lat0,lon0,lat1,lon1))
}


## French data here http://www.galichon.com/codesgeo/
## Other data http://earth-info.nga.mil/gns/html/cntry_files.html

degToRad <- function(val)
{
  return (val * pi / 180.0)
}

## return postal codes nearest to n Kms from a specific code
codesNearToCode <- function(code,kms)
{

  findCode = myGeoData[,2] == code
  latRef = (myGeoData[findCode,5])[1]
  lonRef = (myGeoData[findCode,6])[1]
  distFromCode = distKmRad(latRef,lonRef, myGeoData[,5], myGeoData[,6] )
  return (sort(unique( myGeoData[distFromCode<kms,2] )))

}

## return postal codes nearest to n Kms from a specific code
cityNearToCode <- function(code,kms)
{
  findCode = myGeoData[,2] == code
  latRef = (myGeoData[findCode,5])[1]
  lonRef = (myGeoData[findCode,6])[1]
  distFromCode = distKmRad(latRef,lonRef, myGeoData[,5], myGeoData[,6] )
  return ( myGeoData[distFromCode<kms,1] )

}

# return code of city
zipCode <- function(city)
{

  allCities = toupper(as.character(myGeoData[,1]))
  allCities = gsub("[ \t-]","",allCities)

  res =  myGeoData[allCities ==gsub("[ \t-]","",toupper(city)) ,2]

  return( res)
}

## return cities that have specific zipcode
cities <- function(code)
{

  return (as.character( myGeoData[myGeoData[,2]==code,1] ))
}
