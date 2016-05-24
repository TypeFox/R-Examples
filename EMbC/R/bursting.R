# The EMbC Package for R
# 
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic Bartumeus <fbartu@ceab.csic.es>
#   
#   EMbC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
# 
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.

# binClstPath burst-wise computation

#------------------------------------------------------------------------------------
# John version
#------------------------------------------------------------------------------------

# calculates destination lat and lon along loxodrome given starting point, distance,
# and constant bearing

loxDest = function(lon1deg, lat1deg,  distMeters, brngdeg){
  
  torad = pi/180
  todeg = 180/pi
  
  d = distMeters/earthR()
  lat1 = lat1deg*torad
  lon1 = lon1deg*torad
  brng = brngdeg*torad
  
  dLat = d*cos(brng)
  lat2 = lat1 + dLat;
  dPhi = log(tan(lat2/2+pi/4)/tan(lat1/2+pi/4))
  q = ifelse((dPhi!=0 && dPhi!='-Inf.'), dLat/dPhi, cos(lat1))  
  dLon = d*sin(brng)/q
  
  if (abs(lat2) > pi/2) lat2 = ifelse(lat2>0, pi-lat2, -pi-lat2)
  
  lon2 = (lon1+dLon+pi)%%(2*pi) - pi
  
  return(c(lon2*todeg, lat2*todeg))
}

# finding loxodromic midpoint
findMidPoint <- function(i, breakDf, bCP, myLines){
  
  newLines = coordinates(myLines)[[i]]
  if(breakDf$distance[i]==0){
    lat1 = coordinates(myLines)[[i]][[1]][1,"lat"]
    lon1 = coordinates(myLines)[[i]][[1]][1,"lon"]
    return(list(cbind(lon1, lat1), newLines))}
  
  mpDist = breakDf$distance[i]/2
  theseSegsDists = bCP@dst[breakDf$first[i]:breakDf$last[i]]
  theseSegsCumDists = cumsum(theseSegsDists)
  mpSeg = min(which(theseSegsCumDists > mpDist))
  d = theseSegsDists[mpSeg] - (theseSegsCumDists[mpSeg]- mpDist)
  b = bCP@hdg[breakDf$first[i]:breakDf$last[i]][mpSeg]
  lat1 = coordinates(myLines)[[i]][[1]][mpSeg,"lat"]
  lon1 = coordinates(myLines)[[i]][[1]][mpSeg,"lon"]
  
  onLastSegment <- (mpSeg==nrow(coordinates(myLines)[[i]][[1]]))
  if (onLastSegment){
    lat2 <- coordinates(myLines)[[i+1]][[1]][1,"lat"]
    lon2 <- coordinates(myLines)[[i+1]][[1]][1,"lon"]}
  else{
    lat2 <- coordinates(myLines)[[i]][[1]][mpSeg+1,"lat"]
    lon2 <- coordinates(myLines)[[i]][[1]][mpSeg+1,"lon"]}
  
  if(all(c(lon1,lat1)==c(lon2,lat2))){
    return(list(cbind(lon1, lat1), newLines))}
  else {
    res = t(as.matrix(loxDest(lat1deg=lat1, lon1deg=lon1, distMeters=d, brngdeg=b*180/pi)))
    
    if(onLastSegment||breakDf$distance[i]==0){
      newLines[[1]] =  rbind(coordinates(myLines)[[i]][[1]][1:mpSeg,], res)
    } else{
      newLines[[1]] =  rbind(coordinates(myLines)[[i]][[1]][1:mpSeg,], res, coordinates(myLines)[[i]][[1]][(mpSeg+1):nrow(coordinates(myLines)[[i]][[1]]),])
    }
    colnames(res) = c("lon", "lat")
    return(list(res, newLines))
  }
}

# burst-wise computation
burstIt <- function(bCP,extObjName){
				
			# creating a dataframe showing first and last points for each segment, along with that segment's annotation color.
			brks <- which(diff(bCP@A[1:(bCP@n-1)])!=0)+1
			breakDf = data.frame("A"=bCP@A[c(1,brks)],"first"=c(1,brks),"last"=c(brks,bCP@n),row.names=as.character(1:(length(brks)+1)))  
#			breakDf = data.frame("A"=bCP@A[brks],"first"=c(1,brks[1:length(brks)-1]),"last"=brks,row.names=as.character(1:length(brks)))  
			# now a spatial lines dataframe
			myLines = SpatialLinesDataFrame(SpatialLines(lapply(1:nrow(breakDf), function(i)
			Lines(Line(bCP@pth[breakDf[i,"first"]:breakDf[i,"last"],2:3]), ID=i)),proj4string=CRS("+proj=longlat")), breakDf)  
			# calculating durations. (Note that we subtract 1 from the last index because this is the endpoint, and its duration applies to the subsequent segment)
			breakDf$duration = unlist(lapply(1:nrow(breakDf), function(i) sum(bCP@spn[breakDf$first[i]:(breakDf$last[i]-1)])))
			# calculating distance. (Note that we subtract 1 from the last index because this is the endpoint, and its distance applies to the subsequent segment)
			breakDf$distance = unlist(lapply(1:nrow(breakDf), function(i) sum(bCP@dst[breakDf$first[i]:(breakDf$last[i]-1)])))
			# creating new line list for displaying results so that we can force lines though loxodromix midpoints
#this is not used!!	newLineList = coordinates(myLines)
			# calculating midpoints.
			# NOTE: WITH LARGE DATASETS, THIS MAY BE A LITTLE SLOW.
			# WE COULD EASILY MAKE AN MCAPPLY OPTION TO SPEED IT UP ON MULTICORE SYSTEMS
			midResults = lapply(1:nrow(breakDf), function(i) findMidPoint(i,breakDf,bCP,myLines))
			midResultsMids = lapply(midResults, function(thisResult) thisResult[[1]])
			breakDf$lon = t(array(unlist(midResultsMids), c(2,nrow(breakDf))))[,1]
			breakDf$lat = t(array(unlist(midResultsMids), c(2,nrow(breakDf))))[,2]
			# making a spatial points dataframe of the midpoints
			midPointSPDF = breakDf
			coordinates(midPointSPDF) = ~lon + lat
			# Making spatial lines from new line list
			myMappingLines = SpatialLinesDataFrame(SpatialLines(lapply(1:length(midResults),function(i){
				Lines(Line(midResults[[i]][[2]]),ID=i)}),proj4string=CRS("+proj=longlat")),breakDf)  

			# assigning the two new burst objects to the binClst that was used as input
			bCP@bursted <- TRUE 
			bCP@tracks <- myMappingLines
			bCP@midPoints <- midPointSPDF
			assign(extObjName,bCP,envir=parent.frame())
			return(bCP)}
