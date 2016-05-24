#This is the alternate regridding algorithm

#org_area <- NULL

#Load dev data

#orgCmip5 <- loadCMIP5(path='sampledata/fx/', experiment='historical',
#                      variable='areacella', model='HadGEM2-ES')
#orgLon <- cmip5obj$lon
#orgLat <- cmip5obj$lat
#check (?make) orginal area file
#if(is.null(org_area)){
#   org_area <- calcGridArea(lon=orgLon, lat=orgLat)
#}

#lon_mtr <- matrix(orgLon, nrow=length(orgLon),ncol=length(orgLat), byrow=TRUE)
#lat_mtr <- matrix(orgLat, nrow=length(orgLon),ncol=length(orgLat))

#dump orginal area and values to data frame

#calculate upper and lower bounds of bands
#...wrapping lower bounds are negative
#...wrapping upper bounds are positive

#make projection grid area

#dump grid area to data frame

#for each projected grid calculate overlap percent of each orginal grid
#and save thsi into a column of the transfer matrix
##...if the lower bound is greater then upper bound (lat or lon) => 0
##...if the upper bound is less then the lower bound (lat or lon) => 0
##...otherwise overlap is product of (org-proj)/org for lat and lon

#Apply the transfer matrix to every time-depth/lev

#return regridded matrix
