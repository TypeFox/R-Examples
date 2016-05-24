### R code from vignette source 'IRISSeismic-intro.Rnw'

###################################################
### code chunk number 1: IRISSeismic-intro.Rnw:93-95
###################################################
library(IRISSeismic)
iris <- new("IrisClient")


###################################################
### code chunk number 2: IRISSeismic-intro.Rnw:112-117
###################################################
starttime <- as.POSIXct("2002-04-20", tz="GMT")
endtime <- as.POSIXct("2002-04-21", tz="GMT")
st <- getDataselect(iris,"US","OXF","","BHZ",starttime,endtime)
length(st@traces)
plotUpDownTimes(st, min_signal=1, min_gap=1)


###################################################
### code chunk number 3: IRISSeismic-intro.Rnw:128-129
###################################################
getGaps(st)


###################################################
### code chunk number 4: IRISSeismic-intro.Rnw:135-138
###################################################
parallelLength(st)
parallelMax(st)
parallelSd(st)


###################################################
### code chunk number 5: IRISSeismic-intro.Rnw:145-147
###################################################
tr <- st@traces[[3]]
tr@stats


###################################################
### code chunk number 6: IRISSeismic-intro.Rnw:152-153
###################################################
plot(tr)


###################################################
### code chunk number 7: IRISSeismic-intro.Rnw:175-176
###################################################
slotNames(st)


###################################################
### code chunk number 8: IRISSeismic-intro.Rnw:196-199
###################################################
class(st@url)
class(st@requestedStarttime)
class(st@traces)


###################################################
### code chunk number 9: IRISSeismic-intro.Rnw:209-210
###################################################
slotNames(st@traces[[1]])


###################################################
### code chunk number 10: IRISSeismic-intro.Rnw:246-253
###################################################
as.POSIXct("2010-02-27", tz="GMT") # good
as.POSIXct("2010-02-27 04:00:00", tz="GMT") # good
as.POSIXct("2010-02-27T04:00:00", tz="GMT",
           format="%Y-%m-%dT%H:%M:%OS") # good

as.POSIXct("2010-02-27") # BAD -- no timezone
as.POSIXct("2010-02-27T04:00:00", tz="GMT") # BAD -- no formatting


###################################################
### code chunk number 11: IRISSeismic-intro.Rnw:268-269
###################################################
help("IRISSeismic",package="IRISSeismic")


###################################################
### code chunk number 12: IRISSeismic-intro.Rnw:297-311
###################################################
starttime <- as.POSIXct("2010-02-27", tz="GMT")
endtime <- as.POSIXct("2010-02-28", tz="GMT")
st <- getDataselect(iris,"IU","ANMO","00","BHZ",starttime,endtime)

start2 <- as.POSIXct("2010-02-27 06:40:00", tz="GMT")
end2 <- as.POSIXct("2010-02-27 07:40:00", tz="GMT")

tr1 <- st@traces[[1]]
tr2 <- slice(tr1, start2, end2)

layout(matrix(seq(2)))        # layout a 2x1 matrix
plot(tr1)                     # top
plot(tr2)                     # bottom
layout(1)                     # restore original layout


###################################################
### code chunk number 13: IRISSeismic-intro.Rnw:343-350
###################################################
starttime <- as.POSIXct("2002-04-20", tz="GMT")
endtime <- as.POSIXct("2002-04-21", tz="GMT")
st <- getDataselect(iris,"US","OXF","","BHZ",starttime,endtime)
tr <- st@traces[[3]]
picker <- STALTA(tr,3,30)
threshold <- quantile(picker,0.99999,na.rm=TRUE)
to <- triggerOnset(tr,picker,threshold)


###################################################
### code chunk number 14: IRISSeismic-intro.Rnw:372-382
###################################################
layout(matrix(seq(3)))        # layout a 3x1 matrix
closeup1 <- eventWindow(tr,picker,threshold,3600)
closeup2 <- eventWindow(tr,picker,threshold,600)
plot(tr)
abline(v=to, col='red', lwd=2)
plot(closeup1)
abline(v=to, col='red', lwd=2)
plot(closeup2)
abline(v=to, col='red', lwd=2)
layout(1)                     # restore original layout


###################################################
### code chunk number 15: IRISSeismic-intro.Rnw:400-404
###################################################
starttime <- as.POSIXct("2010-02-27", tz="GMT")
endtime <- as.POSIXct("2010-02-28", tz="GMT")
availability <- getAvailability(iris,"IU","ANMO","*","B??",starttime,endtime)
availability


###################################################
### code chunk number 16: IRISSeismic-intro.Rnw:453-501
###################################################
# Open a connection to IRIS DMC webservices
iris <- new("IrisClient")

# Two days around the "Nisqually Quake"
starttime <- as.POSIXct("2001-02-27", tz="GMT")
endtime <- starttime + 3600 * 24 * 2

# Find biggest seismic event over these two days -- it's the "Nisqually"
events <- getEvent(iris, starttime, endtime, minmag=5.0)
bigOneIndex <- which(events$magnitude == max(events$magnitude))
bigOne <- events[bigOneIndex,]

# Find US stations that are available within 10 degrees of arc of the 
# event location during the hour after the event
start <- bigOne$time
end <- start + 3600
av <- getAvailability(iris, "US", "", "", "BHZ", start, end,
                      latitude=bigOne$latitude, longitude=bigOne$longitude,
                      minradius=0, maxradius=10)
    
# Get the station the furthest East
minLonIndex <- which(av$longitude == max(av$longitude))
snclE <- av[minLonIndex,]

# Get travel times to this station
traveltimes <- getTraveltime(iris, bigOne$latitude, bigOne$longitude, bigOne$depth,
                             snclE$latitude, snclE$longitude)

# Look at the list                             
traveltimes

# Find the P and S arrival times
pArrival <- start + traveltimes$travelTime[traveltimes$phaseName=="P"]
sArrival <- start + traveltimes$travelTime[traveltimes$phaseName=="S"] 

# Get the BHZ signal for this station
st <- getDataselect(iris,snclE$network,snclE$station,
                    snclE$location,snclE$channel,
                    start,end)

# Check that there is only a single trace
length(st@traces)

# Plot the seismic trace and mark the "P" and "S" arrival times
tr <- st@traces[[1]]
plot(tr, subsampling=1) # need subsmpling=1 to add vertical lines with abline()
abline(v=pArrival, col='red')
abline(v=sArrival, col='blue')


