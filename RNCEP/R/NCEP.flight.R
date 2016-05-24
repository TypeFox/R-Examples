NCEP.flight <- function(beg.loc, end.loc, begin.dt, flow.assist='NCEP.Tailwind', fa.args=list(airspeed=12), path='loxodrome', calibrate.dir=FALSE, calibrate.alt=TRUE, cutoff=0, when2stop=list('latitude','longitude',50), levels2consider=c(850,925), hours=12, evaluation.interval=60, id=1, land.if.bad=FALSE, reanalysis2 = FALSE, query=TRUE){

## Assign the appropriate function to f.asst from the equation specified above ##
f.asst <- eval(parse(text=paste(flow.assist)))

## If not querying data, reset levels2consider ##
if(query == FALSE & length(levels2consider) > 1){
	levels2consider <- 1
	warning('levels2consider was reset to an arbitrary value of 1 because data are not queried from NCEP')
	}

## If beg and end lat doesn't differ, suggest that when2stop not include latitude
if(beg.loc[1] == end.loc[1] & any(when2stop == 'latitude')){
	stop("Cannot include 'latitude' in when2stop if beginning and ending latitude are the same")
	}
if(beg.loc[2] == end.loc[2] & any(when2stop == 'longitude')){
	stop("Cannot include 'longitude' in when2stop if beginning and ending longitude are the same")
	}

#########################################
## Create vectors to store output data ##
for(p in levels2consider){
	eval(parse(text=paste("u",p," <- c()", sep='')))
	eval(parse(text=paste("v",p," <- c()", sep='')))
	}
lat <- beg.loc[1]
lon <- beg.loc[2]
dt <- begin.dt
dist2goal <- deg.dist(long1=lon, lat1=lat, long2=end.loc[2], lat2=end.loc[1])
fa <- c()
best <- c()
angle <- c()


##############################################
## Create a one-sided 'is.between' function ##
is.between <- function(x, a, b) {
if(b > a) { x < b } else
if(a > b) { x > b }
}
##

################################
## Begin the loop
loop <- 1
dist.adjusted <- FALSE
failed <- FALSE
quit <- FALSE
while (quit==FALSE){ 
##

#######################################
## Calculate the preferred direction ##
if(calibrate.dir == TRUE | loop == 1){
	if(path == 'loxodrome' | path == 'NCEP.loxodrome'){
		angle[loop] <- NCEP.loxodrome(lat1=lat[loop],lat2=end.loc[1],lon1=lon[loop],lon2=end.loc[2])
		} else
	if(path == 'great.circle' | path == 'NCEP.great.circle'){
		angle[loop] <- earth.bear(long1=lon[loop], lat1=lat[loop], long2=end.loc[2], lat2=end.loc[1])
		} else 
	if(is.numeric(path)){
		angle[loop] <- path 
		if(calibrate.dir == TRUE && loop == 1){
		warning("A new flight direction cannot be calculated if 'path' is numeric.")
		}
		}
		} else { angle[loop] <- angle[loop-1] }
##

#####################################################################################################################
## Specify whether to consider all original pressure levels or restrict to the level chosen on the first iteration ##
if(calibrate.alt == FALSE && loop == 2){
	levels2consider <- best 
	}
##

#####################################################################
## Query the U and V wind components for the location and timestep ##
if(query == TRUE){
for(p in levels2consider){
	if(p == 'surface'){
		eval(parse(text=paste("u", p, "[loop] <- NCEP.interp(variable='uwnd.sig995',level=p, lat=lat[loop], lon=lon[loop], dt=dt[loop], return.units=FALSE, status.bar=FALSE, reanalysis2=reanalysis2)",sep='')))
		##
		eval(parse(text=paste("v", p, "[loop] <- NCEP.interp(variable='vwnd.sig995',level=p, lat=lat[loop], lon=lon[loop], dt=dt[loop], return.units=FALSE, status.bar=FALSE, reanalysis2=reanalysis2)",sep='')))
		} else
	if(p == 'gaussian'){
		eval(parse(text=paste("u", p, "[loop] <- NCEP.interp(variable='uwnd.10m',level=p, lat=lat[loop], lon=lon[loop], dt=dt[loop], return.units=FALSE, status.bar=FALSE, reanalysis2=reanalysis2)",sep='')))
		##
		eval(parse(text=paste("v", p, "[loop] <- NCEP.interp(variable='vwnd.10m',level=p, lat=lat[loop], lon=lon[loop], dt=dt[loop], return.units=FALSE, status.bar=FALSE, reanalysis2=reanalysis2)",sep='')))
		} else
	if(p %in% c(1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10)){
		eval(parse(text=paste("u", p, "[loop] <- NCEP.interp(variable='uwnd',level=as.numeric(p), lat=lat[loop], lon=lon[loop], dt=dt[loop], return.units=FALSE, status.bar=FALSE, reanalysis2=reanalysis2)",sep='')))
		eval(parse(text=paste("v", p, "[loop] <- NCEP.interp(variable='vwnd',level=as.numeric(p), lat=lat[loop], lon=lon[loop], dt=dt[loop], return.units=FALSE, status.bar=FALSE, reanalysis2=reanalysis2)",sep='')))
		} else stop("levels2consider must be 'surface', 'gaussian', or a numeric value indicating a pressure level in the NCEP dataset.  See NCEP.interp for possible pressure levels.") 
	}
	} else
if(query == FALSE){
	for(p in levels2consider){
		eval(parse(text=paste("u", p, "[loop] <- NA",sep='')))
		##
		eval(parse(text=paste("v", p, "[loop] <- NA",sep='')))
		}
	}
##

######################################################################################################################
## Calculate the flow-assistance at each pressure level, determine the best level, and assign that value to flow-assistance ##
if(is.null(fa.args) & loop == 1){
	warning("No arguments passed to the flow-assistance equation, default values (if available) will be used.")
	} else
if(is.list(fa.args) | is.null(fa.args)){
	fas <- c()
	if(query == TRUE){
	for(j in 1:length(levels2consider)){
		fas[j] <- do.call('f.asst', args=c(list(u=eval(parse(text=paste("u",levels2consider[j],"[loop]",sep='')))), list(v=eval(parse(text=paste("v",levels2consider[j],"[loop]",sep='')))), list(direction=angle[loop]), fa.args))$fa
		}
		} else
if(query == FALSE){
	for(j in 1:length(levels2consider)){
		fas[j] <- do.call('f.asst', args=c(list(direction=angle[loop], lat.x=lat[loop], lon.x=lon[loop], dt.x=dt[loop]), fa.args))$fa
		}
		}		
	} else { stop('fa.args must be of type list.') }
##
best.fa.loc <- ifelse(all(is.na(fas) == TRUE), 1, which(fas == max(fas)))
fa[loop] <- fas[best.fa.loc]
best[loop] <- levels2consider[best.fa.loc]
##

############################################################################################################
## Abort flight night if all flow conditions from the first measurement of the night are below the cutoff ##
if(loop == 1 && is.na(fa[loop])) { warning("At take-off, the equation produced no real solution for all pressure levels considered.")
movement <- data.frame(NA)
break }
if(loop == 1 && !is.na(fa[loop]) && fa[loop] <= cutoff) { warning(paste("No pressure levels have sufficient flow-assistance to initiate takeoff, consider adjusting 'cutoff'. Best flow-assistance was ", round(fa[loop], digits=3),".",sep=''))
movement <- data.frame(NA)
break }
if(is.na(fa[loop])) { warning("Flight was aborted because the equation produced no real solution for all pressure levels considered.")
movement[loop,] <- NA
failed <- TRUE
break }
if(land.if.bad && fa[loop] <= cutoff ) { warning(paste("Flight was aborted because no pressure levels have sufficient flow-assistance to continue flight, consider adjusting 'cutoff'. Best flow-assistance was ", round(fa[loop], digits=3),".",sep=''))
movement[loop,] <- NA
failed <- TRUE
break }
##

#################################################################################
## Calculate the forward and sideways motion of the bird including the flow    ##
## Using the specified flow-assistance strategy and the best fa pressure level ##
if(query == TRUE){
if(loop == 1){
movement <- do.call('f.asst', args=c(list(u=eval(parse(text=paste("u",best[loop],"[loop]",sep='')))), list(v=eval(parse(text=paste("v",best[loop],"[loop]",sep='')))), list(direction=angle[loop]), fa.args))
	} else {
movement[loop,] <- do.call('f.asst', args=c(list(u=eval(parse(text=paste("u",best[loop],"[loop]",sep='')))), list(v=eval(parse(text=paste("v",best[loop],"[loop]",sep='')))), list(direction=angle[loop]), fa.args))
	}
} else
if(query == FALSE){
if(loop == 1){
movement <- do.call('f.asst', args=c(list(direction=angle[loop], lat.x=lat[loop], lon.x=lon[loop], dt.x=dt[loop]), fa.args))
	} else {
movement[loop,] <- do.call('f.asst', args=c(list(direction=angle[loop], lat.x=lat[loop], lon.x=lon[loop], dt.x=dt[loop]), fa.args))
	}
} 
	
#################################################################################################################
## Calculate the distance forward and sideways (in meters) that the bird will travel until the next evaluation ##
move.forward <- as.numeric(movement$forward.move[loop])*60 ## Now in meters per minute 
move.side <- as.numeric(movement$side.move[loop])*60 ## Now in meters per minute 
##

###################################################################################################
## Calculate the new bearing and travel distance (in meters) after incorporating flow conditions ##
new.bearing <- angle[loop] + atan2(move.side,move.forward)*(180/pi) 
new.distance <- sqrt(move.side^2+move.forward^2)
##

###################################################################################################################
## Iteratively move the bird to the next location according to its speed, preferred direction, and interaction with the flow ##
latlon <- matrix(NA, ncol=2, nrow=evaluation.interval+1)
latlon[1,] <- c(lat[loop],lon[loop])
for(i in 2:length(latlon[,1])){
	latlon[i,] <- new.lat.long(long=latlon[i-1,2], lat=latlon[i-1,1], bearing=new.bearing, distance=new.distance/1000)
	}
new.loc <- latlon[i,]


#################################################################################################
## See if the bird has passed its specified latitude or longitude between evaluation intervals ##
between.latlon <- cbind(is.between(latlon[,1], beg.loc[1], end.loc[1]), is.between(latlon[,2], beg.loc[2], end.loc[2]))

###################################################################################
## See if the bird has passed through the goal area between evaluation intervals ##
if(any(lapply(when2stop, FUN=is.numeric) == TRUE)){
outside.of.area <- deg.dist(long1=latlon[,2], lat1=latlon[,1], long2=end.loc[2], lat2=end.loc[1]) > when2stop[[which(lapply(when2stop, FUN=is.numeric) == TRUE)]]
}

############################################################################
## Adjust ending location and time spent travelling accordingly if needed ##
if(any(when2stop == 'latitude') && any(between.latlon[,1] == FALSE) && any(when2stop == 'longitude') && any(between.latlon[,2] == FALSE) && loop != 1) {
	lat.pass <- which(between.latlon[,1] == FALSE)[1]
	lon.pass <- which(between.latlon[,2] == FALSE)[1]
	pass <- min(c(lat.pass, lon.pass)) - 1
	new.loc <- latlon[pass,]
	dist.adjusted <- TRUE
	} else
if(any(when2stop == 'latitude') && any(between.latlon[,1] == FALSE) && loop != 1) {
	pass <- which(between.latlon[,1] == FALSE)[1] - 1
	new.loc <- latlon[pass,]
	dist.adjusted <- TRUE
	 } else
if(any(when2stop == 'longitude') && any(between.latlon[,2] == FALSE) && loop != 1) {
	pass <- which(between.latlon[,2] == FALSE)[1] - 1
	new.loc <- latlon[pass,]
	dist.adjusted <- TRUE
	} else
if(any(lapply(when2stop, FUN=is.numeric) == TRUE) && any(outside.of.area == FALSE)){
	pass <- which(outside.of.area == FALSE)[1]
	new.loc <- latlon[pass,]
	dist.adjusted <- TRUE
	}


###################################################
## Advance the loop number to the next iteration ##
loop <- loop+1


################################
## Determine the new datetime ##
if(dist.adjusted == FALSE){
	dt[loop] <- format(strptime(dt[loop-1], format="%Y-%m-%d %H:%M:%S", tz='UTC') + evaluation.interval*60, "%Y-%m-%d %H:%M:%S")
} else
if(dist.adjusted == TRUE){
	dt[loop] <- format(strptime(dt[loop-1], format="%Y-%m-%d %H:%M:%S", tz='UTC') + (pass*60), "%Y-%m-%d %H:%M:%S")
	}

#############################
## Assign the new location ##
lat[loop] <- new.loc[1]
lon[loop] <- new.loc[2]

##################################################
## Assign the new distance to the end goal (km) ##
dist2goal[loop] <- deg.dist(long1=lon[loop], lat1=lat[loop], long2=end.loc[2], lat2=end.loc[1])

##
######################################
## See if the goal has been reached ##
if(dist.adjusted == TRUE) { 
	quit <- TRUE
	}
if(loop != 1 && as.numeric(difftime(as.POSIXct(dt[loop], format="%Y-%m-%d %H:%M:%S", tz='UTC'), as.POSIXct(begin.dt, format="%Y-%m-%d %H:%M:%S", tz='UTC'), units='hours')) >= hours) { 
	quit <- TRUE 
	}
##

} ## Proceed to the next iteration ##


#######################################
## Put the flow data in a data.frame ##
flow.levels <- c()
if(loop == 1 | failed == TRUE){
for(p in levels2consider){ 
	flow.levels <- paste(flow.levels, paste('u_',p,'=c(u',p, '), v_',p,'=c(v',p,')',sep=''), sep=',') 
	}
} else {
for(p in levels2consider){ 
	flow.levels <- paste(flow.levels, paste('u_',p,'=c(u',p, ',NA), v_',p,'=c(v',p,',NA)',sep=''), sep=',') 
	}
}
## Format the datetime variable to datetime format ##
datetime <- strptime(dt, format="%Y-%m-%d %H:%M:%S", tz="UTC")

############################################################
## Put all important variables in a data.frame for output ##
if(loop == 1){
output <- eval(parse(text=paste("data.frame(id=rep(id, length(dt)), datetime=datetime, lat=lat, lon=lon, dist2goal=dist2goal, best=best, flightangle=angle", flow.levels," )", sep='')))
} else
if(failed == TRUE){
output <- eval(parse(text=paste("data.frame(id=rep(id, length(dt)), datetime=datetime, lat=lat, lon=lon, dist2goal=dist2goal, best=best, flightangle=angle", flow.levels,", movement)", sep='')))
} else {
output <- eval(parse(text=paste("data.frame(id=rep(id, length(dt)), datetime=datetime, lat=lat, lon=lon, dist2goal=dist2goal, best=c(best, NA), flightangle=c(angle,NA)", flow.levels,", rbind(movement, NA))", sep='')))
}

return(output)

}
