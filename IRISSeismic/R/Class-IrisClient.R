##
##    S4 classes for communicating with IRIS DMC web services.
##
##    Copyright (C) 2012  Mazama Science, Inc.
##    by Jonathan Callahan, jonathan@mazamascience.com
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

################################################################################
# R classes for an IRIS web servivces request client.
#
# Inspiration for the functionality in this code comes from obspy.iris.client.Client 
# 
#   http://docs.obspy.org/packages/autogen/obspy.iris.client.Client.html
#
#
# Documentation for FDSN Web Services is available here:
#
#   http://www.fdsn.org/webservices/
#
################################################################################

################################################################################
# Class IrisClient
#
# IRIS Web service request client.
#
# Documentation verbatim from obspy.iris.client.Client:
#
#    :type site: str, optional
#    :param site: Base URL of the IRIS Web service (default
#        is ``'http://service.iris.edu/'``).
#    :type user: str, optional
#    :param user: The user name used for authentication with the Web
#        service (default an empty string).
#    :type password: str, optional
#    :param password: A password used for authentication with the Web
#        service (default is an empty string).
#    :type timeout: int, optional
#    :param timeout: Seconds before a connection timeout is raised (default
#        is ``10`` seconds). This works only for Python >= 2.6.x.
#    :type debug: bool, optional
#    :param debug: Enables verbose output (default is ``False``).
#    :type user_agent: str, optional
#    :param user_agent: Sets an client identification string which may be
#        used on server side for statistical analysis (default contains the
#        current module version and basic information about the used
#        operation system, e.g.
#        ``'ObsPy 0.4.7.dev-r2432 (Windows-7-6.1.7601-SP1, Python 2.7.1)'``.
#
################################################################################

#setClass("IrisClient", 
#         # typed slots (aka attributes) for class IrisClient
#         representation(site = "character",
#                        user = "character",
#                        password = "character",
#                        timeout = "integer",
#                        debug = "logical",
#                        useragent = "character"),
#         # default values for slots
#         prototype(site = "http://service.iris.edu",
#                   user = "",
#                   password = "",
#                   timeout = as.integer(10),
#                   debug = FALSE,
#                   useragent = paste("R package IRISSeismic 0.1-3,", R.version.string, ",", version$platform)
#         )
#)

# For now we'll just include the basics:

# REC - check for a user R profile for the IRIS Client site URL
# otherwise, use the default URL
irisSite <- "http://service.iris.edu"
if (Sys.getenv("IrisClient_site") != "") {
	irisSite <- Sys.getenv("IrisClient_site")
}
# check R profile for specific reference to .netrc authentication file
# which allows us to access restricted data sets
# e.g. irisNetrc <- "/home/mustang/mustang.netrc"
# we will use no default here
irisNetrc <- NULL
if (Sys.getenv("IrisClient_netrc") != "") {
	irisNetrc <- Sys.getenv("IrisClient_netrc")
}

# we will want the client using us as a library to be able to identify itself.
# it's easy enough to make this an environment setting as well.
# we will generate a default user agent to represent the metrics calculators.
irisUserAgent <- "unidentified R script"
if (Sys.getenv("IrisClient_agent") != "") {
	irisUserAgent <- Sys.getenv("IrisClient_agent")
}

setClass("IrisClient", 
         # typed slots (aka attributes) for class IrisClient
         representation(site = "character",
                        debug = "logical",
                        useragent = "character"),
         prototype(site = irisSite,
                   debug = FALSE,
                   useragent = paste0("IRISSeismic ",
                                      ifelse("IRISSeismic" %in% rownames(installed.packages()),
                                             installed.packages()["IRISSeismic","Version"], "local code"),
                                      ", RCurl ",
                                      ifelse("RCurl" %in% rownames(installed.packages()),
                                             installed.packages()["RCurl","Version"], "local code"),
                                      ", ", R.version.string,
                                      ", OS ", version$platform,
                                      ", (",irisUserAgent,")"))
)


################################################################################
# getDataselect method returns a Stream object
#
# This implementation matches some of the functionality in the obspy dataselect
# function.
# 
#   http://docs.obspy.org/_modules/obspy/iris/client.html#Client.dataselect
#
################################################################################

if (!isGeneric("getDataselect")) {
  setGeneric("getDataselect", function(obj, network, station, location, channel, starttime, endtime, quality, ignoreEpoch) {
    standardGeneric("getDataselect")
  })
}

getDataselect.IrisClient <- function(obj, network, station, location, channel, starttime, endtime, quality, ignoreEpoch) {

  # allow authenticated access.
  # if we have a netrc definition, then use queryauth to access data
  # else, use standard query call
  url <- NULL
  if (! is.null(irisNetrc)) {
    	url <- paste(obj@site,"/fdsnws/dataselect/1/queryauth?",sep="")
  } else {
    	url <- paste(obj@site,"/fdsnws/dataselect/1/query?",sep="")
  }
  url <- paste(url,"net=", network,sep="")
  url <- paste(url,"&sta=", station,sep="")
  # NOTE:  Locations with blanks must be converted into "--" when creating the URL
  # NOTE:  For getDataselect only, convert "" to "--"
  location <- ifelse(location=="","--",location)
  url <- paste(url,"&loc=", stringr::str_replace(location,"  ","--"),sep="")
  url <- paste(url,"&cha=", channel,sep="")
  url <- paste(url,"&start=", format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&end=", format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&quality=", ifelse(missing(quality),"B",quality),sep="")
  
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # Make authenticated request using a netrc file
  if (! is.null(irisNetrc)) {
      result <- try( dataselectResponse <- RCurl::getBinaryURL(url, useragent=obj@useragent, 
				  netrc=1, netrc.file=irisNetrc),
		  silent=TRUE)
	} else {
  		result <- try( dataselectResponse <- RCurl::getURLContent(url, useragent=obj@useragent),
                 silent=TRUE)
	}
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Bad Request")) {
      stop(paste("getDataselect.IrisClient: Bad Request:",url))
    } else if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getDataselect.IrisClient: Not Found:",url))
    } else {
      stop(paste("getDataselect.IrisClient:",err_msg))
    } 
    
  }  
  
  # No errors so proceed
  
  # Channel metadata is required to properly apply InstrumentSensitivity corrections
  result <- try( channels <- getChannel(obj, network, station, location, channel, starttime, endtime),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Bad Request")) {
      stop(paste("getDataselect.IrisClient: Bad Request:",url))
    } else {
      stop(paste("getDataselect.IrisClient:",err_msg))
    } 
    
  }  
  
  # NOTE:  Sometimes, the station webservice will return multiple records for the same SNCL
  # NOTE:  each with a different scale or starttime.  This still represents a single SNCL.
  # NOTE:  What to do about multiple scales in the following getChannel request?
  # Solution 1: if ignoreEpoch==TRUE, then just take the first epoch presented
  # http://service.iris.edu/fdsnws/station/1/query?net=H2&sta=H2O&loc=00&cha=LHZ&starttime=2001-02-28T18:29:44&endtime=2001-02-28T19:29:44&includerestricted=false&output=text&level=channel
  sncls <- paste(channels$network,channels$station,channels$location,channels$channel)
  if (nrow(channels) > 1 && !ignoreEpoch) {
    stop(paste("getDataselect.IrisClient: Multiple epochs: getChannel returned",length(sncls),"records"))
  } else {
    channelInfo <- channels[1,]
  }
  
  # No errors so proceed

  stream <- miniseed2Stream(dataselectResponse,url,starttime,endtime,
                            channelInfo$instrument,channelInfo$scale,channelInfo$scaleunits)
  
  return(stream)
}

# All arguments specified
setMethod("getDataselect", signature(obj="IrisClient", 
                                     network="character", station="character", location="character",
                                     channel="character", starttime="POSIXct", endtime="POSIXct",
                                     quality="character", ignoreEpoch="logical"), 
          function(obj, network, station, location, channel, starttime, endtime, quality, ignoreEpoch)
            getDataselect.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                     quality, ignoreEpoch))
# quality="missing", use "B" for "Best", ignoreEpoch="missing", default to FALSE
setMethod("getDataselect", signature(obj="IrisClient", 
                                     network="character", station="character", location="character", 
                                     channel="character", starttime="POSIXct", endtime="POSIXct",
                                     quality="missing", ignoreEpoch="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, quality, ignoreEpoch) 
            getDataselect.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                     quality="B", ignoreEpoch=FALSE))
# ignoreEpoch="missing", default to FALSE
setMethod("getDataselect", signature(obj="IrisClient", 
                                     network="character", station="character", location="character", 
                                     channel="character", starttime="POSIXct", endtime="POSIXct",
                                     quality="character", ignoreEpoch="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, quality, ignoreEpoch) 
            getDataselect.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                     quality, ignoreEpoch=FALSE))
# quality="missing", use "B" for "Best"
setMethod("getDataselect", signature(obj="IrisClient", 
                                     network="character", station="character", location="character", 
                                     channel="character", starttime="POSIXct", endtime="POSIXct",
                                     quality="missing", ignoreEpoch="logical"), 
          function(obj, network, station, location, channel, starttime, endtime, quality, ignoreEpoch) 
            getDataselect.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                     quality="B", ignoreEpoch))


################################################################################
# getSNCL method is a convenience wrapper for getDataselect
################################################################################

if (!isGeneric("getSNCL")) {
  setGeneric("getSNCL", function(obj, sncl, starttime, endtime, quality) {
    standardGeneric("getSNCL")
  })
}

getSNCL.IrisClient <- function(obj, sncl, starttime, endtime, quality) {
  parts <- unlist(stringr::str_split(sncl,'\\.'))
  return( getDataselect.IrisClient(obj, parts[1], parts[2], parts[3], parts[4], starttime, endtime, quality) )
}


# All arguments specified
setMethod("getSNCL", signature(obj="IrisClient",
                               sncl="character", starttime="POSIXct", endtime="POSIXct",
                               quality="missing"),
          function(obj, sncl, starttime, endtime, quality)
            getSNCL.IrisClient(obj, sncl, starttime, endtime, quality))
# quality="missing", use "B" for "Best"
setMethod("getSNCL", signature(obj="IrisClient",
                               sncl="character", starttime="POSIXct", endtime="POSIXct",
                               quality="missing"),
          function(obj, sncl, starttime, endtime, quality)
            getSNCL.IrisClient(obj, sncl, starttime, endtime, quality="B"))


################################################################################
# getRotation method returns a list of three Stream objects
#
# Data are obtained from the rotation web service:
# 
#   http://service.iris.edu/irisws/rotation/1/
#
################################################################################

if (!isGeneric("getRotation")) {
  setGeneric("getRotation", function(obj, network, station, location, channelSet, starttime, endtime, processing) {
    standardGeneric("getRotation")
  })
}

getRotation.IrisClient <- function(obj, network, station, location, channelSet, starttime, endtime, processing) {
  
  url <- paste(obj@site,"/irisws/rotation/1/query?",sep="")
  url <- paste(url,"net=", network,sep="")
  url <- paste(url,"&sta=", station,sep="")
  # NOTE:  Locations with blanks must be converted into "--" when creating the URL
  # NOTE:  For getRotation only, convert "" to "--"
  location <- ifelse(location=="","--",location)
  url <- paste(url,"&loc=", stringr::str_replace(location,"  ","--"),sep="")
  url <- paste(url,"&chaset=", channelSet,sep="")
  url <- paste(url,"&start=", format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&end=", format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  
  url <- paste(url,processing,sep="")    
  
  url <- paste(url,"&output=miniseed",sep="")
  
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # Download the rotation serivce response (zip file) to a tempoary file
  temp <- tempfile()
  utils::download.file(url,temp,quiet=TRUE)
  
  # Create a dataframe of zip file contents and discard the metadata file
  zipContents <- utils::unzip(temp,list=TRUE,junkpaths=TRUE)
  metadata_mask <- stringr::str_detect(zipContents$Name,"metadata")
  zipContents <- zipContents[!metadata_mask,]
  
  # Sanity check
  if (nrow(zipContents) != 3) {
    stop(paste("IrisClient.getRotation: rotation service return has",nrow(zipContents),"traces -- 3 expected"))
  }
  
  # Extract the three miniseed records, convert them to Stream objects and store them in streamList
  streamList <- list()
  for (i in seq(3)) {
    con <- unz(temp,zipContents$Name[i],open="rb")
    rawMiniseed <- readBin(con,"raw",zipContents$Length[i])
    
    streamList[[i]] <- miniseed2Stream(rawMiniseed,
                                       url=paste("rotation:",zipContents$Name[i]),
                                       requestedStarttime=starttime,
                                       requestedEndtime=endtime,
                                       sensor="rotation web service")
    #                                      scale=as.integer(NA),
    #                                      scaleunits="") {    
    close(con)
  }
  
  # Remove the temporary file
  unlink(temp)
  
  return(streamList)    
}

# All arguments specified
setMethod("getRotation", signature(obj="IrisClient", 
                                     network="character", station="character", location="character",
                                     channelSet="character", starttime="POSIXct", endtime="POSIXct",
                                     processing="character"), 
          function(obj, network, station, location, channelSet, starttime, endtime, processing) 
            getRotation.IrisClient(obj, network, station, location, channelSet, starttime, endtime, processing))


################################################################################
# getNetwork method returns a dataframe with information from the output
# of the fdsn station web service with "format=text&level=network".
#
#   http://service.iris.edu/fdsnws/station/1/
#
# Example output:
#
# #Network | Description | StartTime | EndTime | TotalStations
# IU|Global Seismograph Network (GSN - IRIS/USGS)|1988-01-01T00:00:00|2500-12-12T23:59:59|254
#
################################################################################

if (!isGeneric("getNetwork")) {
  setGeneric("getNetwork", function(obj, network, station, location, channel,
                                         starttime, endtime, includerestricted,
                                         latitude, longitude, minradius, maxradius) {
    standardGeneric("getNetwork")
  })
}

getNetwork.IrisClient <- function(obj, network, station, location, channel,
                                       starttime, endtime, includerestricted,
                                       latitude, longitude, minradius, maxradius) {
  
  # Parameters common to all 'station' webservice requests
  url <- paste(obj@site,"/fdsnws/station/1/query?",sep="")
  url <- paste(url,"net=",ifelse(network=="","*",network),sep="")
  url <- paste(url,"&sta=",ifelse(station=="","*",station),sep="")
  # NOTE:  Blank locations containing two spaces must be converted to "--"
  location <- ifelse(location=="","*",location)
  url <- paste(url,"&loc=",stringr::str_replace(location,"  ","--"),sep="")
  url <- paste(url,"&cha=",ifelse(channel=="","*",channel),sep="")
  url <- paste(url,"&starttime=",format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&endtime=",format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&includerestricted=",ifelse(includerestricted,"true","false"),sep="")
  url <- paste(url,"&output=text",sep="")
  
  # Add optional geographic constraints if they are passed in
  if (!missing(latitude)) {
    url <- paste(url,"&latitude=",latitude,sep="")      
  }
  if (!missing(longitude)) {
    url <- paste(url,"&longitude=",longitude,sep="")      
  }
  if (!missing(minradius)) {
    url <- paste(url,"&minradius=",minradius,sep="")      
  }
  if (!missing(maxradius)) {
    url <- paste(url,"&maxradius=",maxradius,sep="")      
  }
  
  # Parameters specific to the getNetwork() method
  url <- paste(url,"&level=network",sep="")
  
  # Write debug output
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # Set up the colnames we wish to have in our dataframe
  # #Network | Description | StartTime | EndTime | TotalStations
  colNames <- c("network","description","starttime","endtime","totalstations")
  colClasses <- c(rep("character",4),"numeric")
  
  # Make webservice request
  # NOTE:  Be sure to set na.strings="" as "NA" is a valid network name  
  #result <- try( DF <- utils::read.delim(url,sep="|",col.names=colNames,colClasses=colClasses,na.strings=""),
  #               silent=TRUE)
  # REC -- going with RCurl variation for proper file closure
  result <- try( DF <- utils::read.delim(textConnection(RCurl::getURL(url,useragent=obj@useragent)),sep="|",col.names=colNames,colClasses=colClasses,na.strings=""),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"204 No Content")) {
      stop(paste("getNetwork.IrisClient: No Content:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getNetwork.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getNetwork.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  # Convert time strings
  DF$starttime <- as.POSIXct(DF$starttime, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  DF$endtime <- as.POSIXct(DF$endtime, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  
  # Return a dataframe with rows ordered by network
  return(DF[order(DF$network),])
}

# All arguments specified
setMethod("getNetwork", signature(obj="IrisClient", 
                                       network="character", station="character", location="character",
                                       channel="character", starttime="POSIXct", endtime="POSIXct",
                                       includerestricted="logical",
                                       latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getNetwork.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                  includerestricted, latitude, longitude, minradius, maxradius))

# includerestricted="missing", use FALSE
setMethod("getNetwork", signature(obj="IrisClient", 
                                       network="character", station="character", location="character",
                                       channel="character", starttime="POSIXct", endtime="POSIXct",
                                       includerestricted="missing",
                                       latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getNetwork.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                  FALSE, latitude, longitude, minradius, maxradius))


################################################################################
# getStation method returns a dataframe with information from the output
# of the fdsn station web service with "format=text&level=station".
#
#   http://service.iris.edu/fdsnws/station/1/
#
# Example output:
#
# #Network | Station | Latitude | Longitude | Elevation | SiteName | StartTime | EndTime 
# IU|COR|44.5855|-123.3046|110.0|Corvallis, Oregon, USA|2009-09-26T00:00:00|2599-12-31T23:59:59
#
################################################################################

if (!isGeneric("getStation")) {
  setGeneric("getStation", function(obj, network, station, location, channel,
                                    starttime, endtime, includerestricted,
                                    latitude, longitude, minradius, maxradius) {
    standardGeneric("getStation")
  })
}

getStation.IrisClient <- function(obj, network, station, location, channel,
                                  starttime, endtime, includerestricted,
                                  latitude, longitude, minradius, maxradius) {
  
  # Parameters common to all 'station' webservice requests
  url <- paste(obj@site,"/fdsnws/station/1/query?",sep="")
  url <- paste(url,"net=",ifelse(network=="","*",network),sep="")
  url <- paste(url,"&sta=",ifelse(station=="","*",station),sep="")
  # NOTE:  Blank locations containing two spaces must be converted to "--"
  location <- ifelse(location=="","*",location)
  url <- paste(url,"&loc=",stringr::str_replace(location,"  ","--"),sep="")
  url <- paste(url,"&cha=",ifelse(channel=="","*",channel),sep="")
  url <- paste(url,"&starttime=",format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&endtime=",format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&includerestricted=",ifelse(includerestricted,"true","false"),sep="")
  url <- paste(url,"&output=text",sep="")
  
  # Add optional geographic constraints if they are passed in
  if (!missing(latitude)) {
    url <- paste(url,"&latitude=",latitude,sep="")      
  }
  if (!missing(longitude)) {
    url <- paste(url,"&longitude=",longitude,sep="")      
  }
  if (!missing(minradius)) {
    url <- paste(url,"&minradius=",minradius,sep="")      
  }
  if (!missing(maxradius)) {
    url <- paste(url,"&maxradius=",maxradius,sep="")      
  }
  
  # Parameters specific to the getStation() method
  url <- paste(url,"&level=station",sep="")
  
  # Write debug output
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # Set up the colnames we wish to have in our dataframe
  # Network | Station | Latitude | Longitude | Elevation | SiteName | StartTime | EndTime 
  colNames <- c("network","station","latitude","longitude","elevation","sitename","starttime","endtime")
  colClasses <- c(rep("character",2),rep("numeric",3),rep("character",3))
  
  # Make webservice request
  # NOTE:  Be sure to set na.strings="" as "NA" is a valid network name  
  #result <- try( DF <- utils::read.delim(url,sep="|",col.names=colNames,colClasses=colClasses,na.strings=""),
  #               silent=TRUE)
  # REC - going with RCurl variation to ensure proper file closure
  result <- try( DF <- utils::read.delim(textConnection(RCurl::getURL(url,useragent=obj@useragent)),sep="|",col.names=colNames,colClasses=colClasses,na.strings=""),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"204 No Content")) {
      stop(paste("getStation.IrisClient: No Content:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getStation.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getStation.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  # Convert time strings
  DF$starttime <- as.POSIXct(DF$starttime, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  DF$endtime <- as.POSIXct(DF$endtime, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  
  # Return dataframe with rows ordered by network.station
  netsta <- paste(DF$network,DF$station,sep='.')
  return(DF[order(netsta),])
}

# All arguments specified
setMethod("getStation", signature(obj="IrisClient", 
                                  network="character", station="character", location="character",
                                  channel="character", starttime="POSIXct", endtime="POSIXct",
                                  includerestricted="logical",
                                  latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getStation.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                  includerestricted, latitude, longitude, minradius, maxradius))

# includerestricted="missing", use FALSE
setMethod("getStation", signature(obj="IrisClient", 
                                  network="character", station="character", location="character",
                                  channel="character", starttime="POSIXct", endtime="POSIXct",
                                  includerestricted="missing",
                                  latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getStation.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                  FALSE, latitude, longitude, minradius, maxradius))


################################################################################
# getChannel method returns a dataframe with information from the output
# of the fdsn station web service with "format=text&level=channel".
#
#   http://service.iris.edu/fdsnws/station/1/
#
# Example output:
#
# #Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip | Instrument | Scale | ScaleFreq | ScaleUnits | SampleRate | StartTime | EndTime
# IU|COR|10|LHZ|44.5855|-123.3046|110.0|0.0|0.0|-90.0|Trillium 240 broad band|1.98775E9|0.02|M/S|1.0|2011-01-15T07:00:00|2012-03-08T08:43:00
#
################################################################################

if (!isGeneric("getChannel")) {
  setGeneric("getChannel", function(obj, network, station, location, channel,
                                    starttime, endtime, includerestricted,
                                    latitude, longitude, minradius, maxradius) {
    standardGeneric("getChannel")
  })
}

getChannel.IrisClient <- function(obj, network, station, location, channel,
                                  starttime, endtime, includerestricted,
                                  latitude, longitude, minradius, maxradius) {
  
  # Parameters common to all 'station' webservice requests
  url <- paste(obj@site,"/fdsnws/station/1/query?",sep="")
  url <- paste(url,"net=",ifelse(network=="","*",network),sep="")
  url <- paste(url,"&sta=",ifelse(station=="","*",station),sep="")
  # NOTE:  Blank locations containing two spaces must be converted to "--"
  location <- ifelse(location=="","*",location)
  url <- paste(url,"&loc=",stringr::str_replace(location,"  ","--"),sep="")
  url <- paste(url,"&cha=",ifelse(channel=="","*",channel),sep="")
  url <- paste(url,"&starttime=",format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&endtime=",format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&includerestricted=",ifelse(includerestricted,"true","false"),sep="")
  url <- paste(url,"&output=text",sep="")
  
  # Add optional geographic constraints if they are passed in
  if (!missing(latitude)) {
    url <- paste(url,"&latitude=",latitude,sep="")      
  }
  if (!missing(longitude)) {
    url <- paste(url,"&longitude=",longitude,sep="")      
  }
  if (!missing(minradius)) {
    url <- paste(url,"&minradius=",minradius,sep="")      
  }
  if (!missing(maxradius)) {
    url <- paste(url,"&maxradius=",maxradius,sep="")      
  }
  
  # Parameters specific to the getChannel() method
  url <- paste(url,"&level=channel",sep="")
  
  # Write debug output
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # Set up the colnames we wish to have in our dataframe
  # Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip |
  # Instrument | Scale | ScaleFreq | ScaleUnits | SampleRate | StartTime | EndTime
  colNames <- c("network","station","location","channel","latitude","longitude","elevation","depth","azimuth","dip",
                "instrument","scale","scalefreq","scaleunits","samplerate","starttime","endtime")
  colClasses <- c(rep("character",4),rep("numeric",6),"character",rep("numeric",2),"character","numeric",rep("character",2))
  
  # Make webservice request
  # NOTE:  Be sure to set na.strings="" as "NA" is a valid network name  
  #result <- try( DF <- utils::read.delim(url,sep="|",col.names=colNames,colClasses=colClasses,na.strings=""),
  #               silent=TRUE)
  # REC -- going with RCurl variation to ensure proper file closure
  result <- try( DF <- utils::read.delim(textConnection(RCurl::getURL(url,useragent=obj@useragent)),sep="|",col.names=colNames,colClasses=colClasses,na.strings=""),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {

    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"204 No Content")) {
      stop(paste("getChannel.IrisClient: No Content:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getChannel.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getChannel.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  # Convert "  " location codes back into the "" that is used in the miniSEED record
  DF$location <- stringr::str_replace(DF$location,"  ","")
  # new station web service now returns blank loc codes, which convert to NA
  # switch NA's to "" as well
  DF$location[is.na(DF$location)] <- ""
  
  # Convert time strings
  DF$starttime <- as.POSIXct(DF$starttime, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  DF$endtime <- as.POSIXct(DF$endtime, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  
  # Add a snclId column
  DF$snclId <- paste(DF$network,DF$station,DF$location,DF$channel,sep=".")
  
  # Return dataframe with rows ordered by snclId
  return(DF[order(DF$snclId),])
}

# All arguments specified
setMethod("getChannel", signature(obj="IrisClient", 
                                  network="character", station="character", location="character",
                                  channel="character", starttime="POSIXct", endtime="POSIXct",
                                  includerestricted="logical",
                                  latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getChannel.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                  includerestricted, latitude, longitude, minradius, maxradius))

# includerestricted="missing", use FALSE
setMethod("getChannel", signature(obj="IrisClient", 
                                  network="character", station="character", location="character",
                                  channel="character", starttime="POSIXct", endtime="POSIXct",
                                  includerestricted="missing",
                                  latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getChannel.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                  FALSE, latitude, longitude, minradius, maxradius))


################################################################################
# getAvailability method returns a dataframe with information from the output
# of the fdsn station web service with "format=text&level=channel".
# With additional parameters, this webservice returns information on all
# matching SNCLs that have available data.
#
# The fdsnws/station/availability web service will return space characters for location
# codes that are SPACE SPACE.
#
#   http://service.iris.edu/fdsnws/station/1/
#
# #Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip | Instrument | Scale | ScaleFreq | ScaleUnits | SampleRate | StartTime | EndTime
# CU|ANWB|00|LHZ|17.66853|-61.78557|39.0|0.0|0.0|-90.0|Streckeisen STS-2 Standard-gain|2.43609E9|0.05|M/S|1.0|2010-02-10T18:35:00|2599-12-31T23:59:59
#
################################################################################

if (!isGeneric("getAvailability")) {
  setGeneric("getAvailability", function(obj, network, station, location, channel,
                                         starttime, endtime, includerestricted,
                                         latitude, longitude, minradius, maxradius) {
    standardGeneric("getAvailability")
  })
}

getAvailability.IrisClient <- function(obj, network, station, location, channel,
                                       starttime, endtime, includerestricted,
                                       latitude, longitude, minradius, maxradius) {

  # Parameters common to all 'station' webservice requests
  url <- paste(obj@site,"/fdsnws/station/1/query?",sep="")
  url <- paste(url,"net=",ifelse(network=="","*",network),sep="")
  url <- paste(url,"&sta=",ifelse(station=="","*",station),sep="")
  # NOTE:  Blank locations containing two spaces must be converted to "--"
  location <- ifelse(location=="","*",location)
  url <- paste(url,"&loc=",stringr::str_replace(location,"  ","--"),sep="")
  url <- paste(url,"&cha=",ifelse(channel=="","*",channel),sep="")
  url <- paste(url,"&starttime=",format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&endtime=",format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&includerestricted=",ifelse(includerestricted,"true","false"),sep="")
  url <- paste(url,"&output=text",sep="")
  
  # Add optional geographic constraints if they are passed in
  if (!missing(latitude)) {
    url <- paste(url,"&latitude=",latitude,sep="")      
  }
  if (!missing(longitude)) {
    url <- paste(url,"&longitude=",longitude,sep="")      
  }
  if (!missing(minradius)) {
    url <- paste(url,"&minradius=",minradius,sep="")      
  }
  if (!missing(maxradius)) {
    url <- paste(url,"&maxradius=",maxradius,sep="")      
  }
    
  # Parameters specific to the getAvailability() method
  url <- paste(url,"&includeavailability=true",sep="")
  url <- paste(url,"&matchtimeseries=true",sep="")
  url <- paste(url,"&level=channel",sep="")

  # Write debug output
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # Set up the colnames we wish to have in our dataframe
  # Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip |
  # Instrument | Scale | ScaleFreq | ScaleUnits | SampleRate | StartTime | EndTime
  colNames <- c("network","station","location","channel","latitude","longitude","elevation","depth","azimuth","dip",
                "instrument","scale","scalefreq","scaleunits","samplerate","starttime","endtime")
  colClasses <- c(rep("character",4),rep("numeric",6),"character",rep("numeric",2),"character","numeric",rep("character",2))
  
  # Make webservice request
  # NOTE:  Be sure to set na.strings="" as "NA" is a valid network name
  #result <- try( DF <- utils::read.delim(url,sep="|",col.names=colNames,colClasses=colClasses,na.strings=""),
  #			silent=TRUE)
  # REC -- going with RCurl variation to ensure proper file closure
  result <- try( DF <- utils::read.delim(textConnection(RCurl::getURL(url,useragent=obj@useragent)),sep="|",col.names=colNames,colClasses=colClasses,na.strings=""),
  			silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"204 No Content")) {
      stop(paste("getAvailability.IrisClient: No Content:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getAvailability.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getAvailability.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed

  # REC Feb 2014
  # Filter out unwanted channels -- first to come to mind is LOG, but others get added by recommendation.
  DF <- DF[!(DF$channel %in% c("LOG", "ACE", "OCF")),]
  
  # Convert "  " location codes back into the "" that is used in the miniSEED record
  DF$location <- stringr::str_replace(DF$location,"  ","")
  # new station web service now returns blank loc codes, which convert to NA
  # switch NA's to "" as well
  DF$location[is.na(DF$location)] <- ""

  # Convert time strings
  DF$starttime <- as.POSIXct(DF$starttime, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  DF$endtime <- as.POSIXct(DF$endtime, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  
  # Add a snclId column
  DF$snclId <- paste(DF$network,DF$station,DF$location,DF$channel,sep=".")
  
  # Return dataframe with rows ordered by snclId
  return(DF[order(DF$snclId),])
}

# All arguments specified
setMethod("getAvailability", signature(obj="IrisClient", 
                                       network="character", station="character", location="character",
                                       channel="character", starttime="POSIXct", endtime="POSIXct",
                                       includerestricted="logical",
                                       latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getAvailability.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                       includerestricted, latitude, longitude, minradius, maxradius))

# includerestricted="missing", use FALSE
setMethod("getAvailability", signature(obj="IrisClient", 
                                       network="character", station="character", location="character",
                                       channel="character", starttime="POSIXct", endtime="POSIXct",
                                       includerestricted="missing",
                                       latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getAvailability.IrisClient(obj, network, station, location, channel, starttime, endtime, 
                                       FALSE, latitude, longitude, minradius, maxradius))


################################################################################
# getUnavailability method returns a dataframe with those channels that are found
# in the getChannel() dataframe but not in the getAvailability() dataframe.
################################################################################

if (!isGeneric("getUnavailability")) {
  setGeneric("getUnavailability", function(obj, network, station, location, channel,
                                           starttime, endtime, includerestricted,
                                           latitude, longitude, minradius, maxradius) {
    standardGeneric("getUnavailability")
  })
}

getUnavailability.IrisClient <- function(obj, network, station, location, channel,
                                         starttime, endtime, includerestricted,
                                         latitude, longitude, minradius, maxradius) {
  
  # Get all channels
  c <- getChannel(obj, network, station, location, channel,
                  starttime, endtime, includerestricted,
                  latitude, longitude, minradius, maxradius)
  
  # Get available channels
  a <- getAvailability(obj, network, station, location, channel,
                       starttime, endtime, includerestricted,
                       latitude, longitude, minradius, maxradius)
  
  # Create unique record identifiers (SNCL is enough)
  c_sncls <- paste(c$network,c$station,c$location,c$channel,sep=".")
  a_sncls <- paste(a$network,a$station,a$location,a$channel,sep=".")
  
  unavailable_mask <- !(c_sncls %in% a_sncls)
  return(c[unavailable_mask,])
}
  
# All arguments specified
setMethod("getUnavailability", signature(obj="IrisClient", 
                                         network="character", station="character", location="character",
                                         channel="character", starttime="POSIXct", endtime="POSIXct",
                                         includerestricted="logical",
                                         latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getUnavailability.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                         includerestricted, latitude, longitude, minradius, maxradius))

# includerestricted="missing", use FALSE
setMethod("getUnavailability", signature(obj="IrisClient", 
                                         network="character", station="character", location="character",
                                         channel="character", starttime="POSIXct", endtime="POSIXct",
                                         includerestricted="missing",
                                         latitude="ANY", longitude="ANY", minradius="ANY", maxradius="ANY"), 
          function(obj, network, station, location, channel, starttime, endtime,
                   includerestricted, latitude, longitude, minradius, maxradius) 
            getUnavailability.IrisClient(obj, network, station, location, channel, starttime, endtime, 
                                         FALSE, latitude, longitude, minradius, maxradius))


################################################################################
# getEvalresp method returns instrument response data from the evalresp webservice:
#
#   http://service/irisws/evalresp/1/
#
################################################################################

if (!isGeneric("getEvalresp")) {
  setGeneric("getEvalresp", function(obj, network, station, location, channel, time,
                                     minfreq, maxfreq, nfreq, units, output) {
    standardGeneric("getEvalresp")
  })
}

getEvalresp.IrisClient <- function(obj, network, station, location, channel, time, 
                                   minfreq, maxfreq, nfreq, units, output) {
  
  # TODO:  getEvalresp should test for wildcards in network, station, location, channel
  # TODO:  and return an error if any are found.
  url <- paste(obj@site,"/irisws/evalresp/1/query?",sep="")
  url <- paste(url,"net=",network,sep="")
  url <- paste(url,"&sta=",station,sep="")
  # NOTE:  Locations with blanks must be converted into "--" when creating the URL
  # NOTE:  For getEvalresp, convert "" to "--"
  location <- ifelse(location=="","--",location)
  url <- paste(url,"&loc=", stringr::str_replace(location,"  ","--"),sep="")
  url <- paste(url,"&cha=",channel,sep="")
  # NOTE:  a single 'time' parameter is used rather than 'starttime' and 'endtime'
  url <- paste(url,"&time=",format(time,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="") # ws_evalresp requires "T" format
  if (missing(output)) {
    output <- "fap"
  }
  url <- paste(url,"&output=",output,sep="")      
  
  # Add optional arguments if they are found
  if (!missing(minfreq)) {
    url <- paste(url,"&minfreq=",minfreq,sep="")      
  }
  if (!missing(maxfreq)) {
    url <- paste(url,"&maxfreq=",maxfreq,sep="")      
  }
  if (!missing(nfreq)) {
    url <- paste(url,"&nfreq=",nfreq,sep="")      
  }
  if (!missing(units)) {
    url <- paste(url,"&units=",units,sep="")      
  }
  
  # Write debug output
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # Handle requests for either freq-amp-phase or freq-real-imag
  if (output == "fap") {
    colNames <- c("freq","amp","phase")
  } else if (output == "cs") {
    colNames <- c("freq","real","imag")
  } else {
    stop(paste("getEvalresp.IrisClient: bad output arg = '",output,"' -- must be 'fap' or 'cs'"))
  }
  
  # Conversion of URL into a data frame is a single line with utils::read.table().
  # REC Jan 2014 -- generate a connection handle for url so we can directly close on error.
  # -- use RCurl variation in read.table call so that failed connections close properly
  #result <- try ( DF <- utils::read.table(con <- url(url), col.names=colNames),
  #                silent=TRUE )
  result <- try ( DF <- utils::read.table(textConnection(RCurl::getURL(url,useragent=obj@useragent)), col.names=colNames),
				  silent=TRUE )
  
  # Handle error response
  if (class(result) == "try-error" ) {
	#
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getEvalresp.IrisClient: URL Not Found:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getEvalresp.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getEvalresp.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  # Return the dataframe, guaranteeing that it is ordered by frequency
    
  return(DF[order(DF$freq),])
}

# All required arguments specified
setMethod("getEvalresp", signature(obj="IrisClient", 
                                   network="character", station="character", location="character",
                                   channel="character", time="POSIXct",
                                   minfreq="ANY", maxfreq="ANY", nfreq="ANY", units="ANY",
                                   output="ANY"), 
          function(obj, network, station, location, channel, time, 
                   minfreq, maxfreq, nfreq, units, output) 
            getEvalresp.IrisClient(obj, network, station, location, channel, time, 
                                   minfreq, maxfreq, nfreq, units, output))


################################################################################
# getEvent method returns seismic event data from the event webservice:
#
#   http://service.iris.edu/fdsnws/event/1/
#
# TODO:  The getEvent method could be fleshed out with a more complete list
# TODO:  of arguments to be used as ws-event parameters.
################################################################################

if (!isGeneric("getEvent")) {
  setGeneric("getEvent", function(obj, starttime, endtime, minmag, maxmag, magtype,
                                  mindepth, maxdepth) {
    standardGeneric("getEvent")
  })
}

getEvent.IrisClient <- function(obj, starttime, endtime, minmag, maxmag, magtype,
                                mindepth, maxdepth) {
  
  url <- paste(obj@site,"/fdsnws/event/1/query?",sep="")
  url <- paste(url,"&starttime=",format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&endtime=",format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="")
  url <- paste(url,"&output=text",sep="")
  
  # Add optional arguments if they are non-null
  if (!missing(minmag)) {
    url <- paste(url,"&minmag=",minmag,sep="")      
  }
  if (!missing(maxmag)) {
    url <- paste(url,"&maxmag=",maxmag,sep="")      
  }
  if (!missing(magtype)) {
    url <- paste(url,"&magtype=",magtype,sep="")      
  }
  if (!missing(mindepth)) {
    url <- paste(url,"&mindepth=",mindepth,sep="")      
  }
  if (!missing(maxdepth)) {
    url <- paste(url,"&maxdepth=",maxdepth,sep="")      
  }
  
  # Write debug output
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # http://service.iris.edu/fdsnws/event/1/query?starttime=2013-02-01T00:00:00&endtime=2013-02-02T00:00:00&minmag=5&output=text
  #
  # #EventID | Time | Latitude | Longitude | Depth | Author | Catalog | Contributor | ContributorID | MagType | Magnitude | MagAuthor | EventLocationName
  # 4075900|2013-02-01T22:18:33|-11.12|165.378|10.0|NEIC|NEIC PDE|NEIC PDE-Q||MW|6.4|GCMT|SANTA CRUZ ISLANDS
  
  # Assign column names and classes for the returned data
  colNames <- c("eventId","time","latitude","longitude","depth","author","cCatalog","contributor",
                "contributorId","magType","magnitude","magAuthor","eventLocationName")
  colClasses <- c("character","character","numeric","numeric","numeric","factor","factor","factor",
                  "factor","factor","numeric","factor","character")
    
  # Conversion of URL into a data frame is a single line with read.table().
  #result <- try ( DF <- utils::read.table(url, sep="|", col.names=colNames, colClasses=colClasses),
  #                silent=TRUE )
  # REC - going with RCurl variation to ensure proper file closure
  result <- try ( DF <- utils::read.table(textConnection(RCurl::getURL(url,useragent=obj@useragent)), sep="|", col.names=colNames, colClasses=colClasses),
                  silent=TRUE )
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg, "cannot open the connection")) {
      stop(paste("getEvent.IrisClient: No Data Found")) 
      # NOTE:  Only "cannot open the connection" has been seen with fdsn web services
      
    } else if (stringr::str_detect(err_msg, "Not Found")) {
      stop(paste("getEvent.IrisClient: URL Not Found:",url))
    } else if (stringr::str_detect(err_msg, "couldn't connect to host")) {
      stop(paste("getEvent.IrisClient: couldn't connect to host"))
    } else if (stringr::str_detect(err_msg, "No Data Found")) {
      stop(paste("getEvent.IrisClient: No Data Found"))
    } else {
      stop(paste("getEvent.IrisClient:",err_msg))
    } 
    
  }
  
  # TODO:  improve error trapping so this last check is not needed
  # Last check to make sure DF is defined
  if (!exists('DF')) {
    stop(paste("getEvent.IrisClient: No Data Found"))            
  }
    
  # No errors so proceed
  
  # Final data conversion
  DF$time <- as.POSIXct(DF$time, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  
  # Return dataframe with rows ordered by time
  return(DF[order(DF$time),])
}
  
# All arguments specified
setMethod("getEvent", signature(obj="IrisClient", 
                                   starttime="POSIXct", endtime="POSIXct",
                                   minmag="ANY", maxmag="ANY", magtype="ANY",
                                   mindepth="ANY", maxdepth="ANY"), 
          function(obj, starttime, endtime, minmag, maxmag, magtype, mindepth, maxdepth) 
            getEvent.IrisClient(obj, starttime, endtime, minmag, maxmag, magtype, mindepth, maxdepth))


################################################################################
# getTraveltime method returns a dataframe with information from the traveltime
# webservice:
#
#   http://service.iris.edu/irisws/traveltime/1/
#
################################################################################

if (!isGeneric("getTraveltime")) {
  setGeneric("getTraveltime", function(obj, latitude, longitude, depth, staLatitude, staLongitude) {
    standardGeneric("getTraveltime")
  })
}

getTraveltime.IrisClient <- function(obj, latitude, longitude, depth, staLatitude, staLongitude) {
  
  # Create URL arguments from incoming parameters
  evloc <- paste("[", latitude, ",", longitude, "]", sep="")
  staloc <- paste("[", staLatitude, ",", staLongitude, "]", sep="")
  
  # Assemble URL
  url <- paste(obj@site,"/irisws/traveltime/1/query?",sep="")
  url <- paste(url,"&evloc=",evloc,sep="")
  url <- paste(url,"&evdepth=",depth,sep="")
  url <- paste(url,"&staloc=", staloc, sep="")
  
  # Write debug output
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # http://service.iris.edu/irisws/traveltime/1/query?evloc=[-11.12,165.378]&evdepth=10.0&staloc=[-30.4183,151.6293]
  #
  #   Model: iasp91
  #   Distance   Depth   Phase   Travel    Ray Param  Takeoff  Incident  Purist    Purist
  #     (deg)     (km)   Name    Time (s)  p (s/deg)   (deg)    (deg)   Distance   Name 
  #   -----------------------------------------------------------------------------------
  #   23.14    10.0   P        306.34    10.542     33.42    33.36    23.14   = P
  #   23.14    10.0   P        306.88     9.162     28.60    28.55    23.14   = P
  #   23.14    10.0   P        308.47     9.679     30.38    30.32    23.14   = P
  #   23.14    10.0   PcP      534.45     2.082      6.25     6.24    23.14   = PcP
  #   23.14    10.0   S        558.79    16.323     29.61    29.55    23.14   = S
  #   23.14    10.0   S        560.09    19.178     35.48    35.42    23.14   = S
  #   23.14    10.0   S        563.00    17.678     32.35    32.29    23.14   = S
  #   23.14    10.0   S        575.53    23.885     46.29    46.20    23.14   = S
  #   23.14    10.0   S        575.59    23.797     46.07    45.98    23.14   = S
  #   23.14    10.0   ScS      978.54     3.844      6.68     6.67    23.14   = ScS
  #   23.14    10.0   PKiKP    998.81     0.511      1.53     1.53    23.14   = PKiKP
  #   23.14    10.0   SKiKS   1422.54     0.570      0.99     0.99    23.14   = SKiKS
  
  # Assign column names and classes for the returned data
  colNames <- c("distance", "depth", "phaseName", "travelTime", "rayParam",
                "takeoff", "incident", "puristDistance", "dummy", "puristName")
  colClasses <- c("numeric", "numeric", "factor", "numeric", "numeric",
                  "numeric", "numeric", "numeric", "character", "character")
  
  # Conversion of URL into a data frame is a single line with read.table().
  #result <- try ( returnValue <- utils::read.table(url, skip=4, col.names=colNames, colClasses=colClasses),
  #                silent=TRUE )
  # REC -- going with RCurl variation for better file closure
  result <- try ( returnValue <- utils::read.table(textConnection(RCurl::getURL(url,useragent=obj@useragent)), skip=4, col.names=colNames, colClasses=colClasses),
                  silent=TRUE )
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getTraveltime.IrisClient: URL Not Found:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getTraveltime.IrisClient: couldn't connect to host"))
    } else if (stringr::str_detect(err_msg,"No Data Found")) {
      stop(paste("getTraveltime.IrisClient: No Data Found"))
    } else {
      stop(paste("getTraveltime.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  # Only use the first seven columns
  DF <- returnValue[,c(1:8,10)]
  
  # Return dataframe with rows ordered by travelTime
  return(DF[order(DF$travelTime),])
}

# All arguments specified
setMethod("getTraveltime", signature(obj="IrisClient", 
                                latitude="numeric", longitude="numeric", depth="numeric",
                                staLatitude="numeric", staLongitude="numeric"), 
          function(obj, latitude, longitude, depth, staLatitude, staLongitude) 
            getTraveltime.IrisClient(obj, latitude, longitude, depth, staLatitude, staLongitude))


################################################################################
# getDistaz method returns a dataframe with information from the distaz
# webservice:
#
#   http://service.iris.edu/irisws/distaz/1/
#
################################################################################

if (!isGeneric("getDistaz")) {
  setGeneric("getDistaz", function(obj, latitude, longitude, staLatitude, staLongitude) {
    standardGeneric("getDistaz")
  })
}

getDistaz.IrisClient <- function(obj, latitude, longitude, staLatitude, staLongitude) {
  
  # Assemble URL
  url <- paste(obj@site,"/irisws/distaz/1/query?",sep="")
  url <- paste(url,"&evtlat=",latitude,sep="")
  url <- paste(url,"&evtlon=",longitude,sep="")
  url <- paste(url,"&stalat=", staLatitude, sep="")
  url <- paste(url,"&stalon=", staLongitude, sep="")
  
  # Write debug output
  if (obj@debug) {
    write(paste("<debug>URL =",url), stdout())
  }
  
  # Get data from distaz web service
  # NOTE:  RCurl::getURLContent returns a binary object based on the "resulting HTTP  header's Content-Type field."
  # NOTE:  Use RCurl::getURL to return data as character.
  result <- try( distazXml <- RCurl::getURL(url, useragent=obj@useragent),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getDistaz.IrisClient: URL Not Found:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getDistaz.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getDistaz.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  # http://service.iris.edu/irisws/distaz/1/query?stalat=0.0&stalon=0.0&evtlat=15.0&evtlon=0.0
  #
  # <DistanceAzimuth>
  #   <azimuth>180.0</azimuth>
  #    <backAzimuth>0.0</backAzimuth>
  #    <distance>14.90407</distance>
  #  </DistanceAzimuth>
  
  # See this example for magic with the XML package:
  #
  #  http://www.omegahat.org/RSXML/gettingStarted.html
  
  doc <- XML::xmlRoot(XML::xmlTreeParse(distazXml)) 
  namedValueStrings <- XML::xmlSApply(doc,XML::xmlValue)
  DF <- as.data.frame(t(as.numeric(namedValueStrings)))
  names(DF) <- names(namedValueStrings)
  
  # Return dataframe with rows ordered by travelTime
  return(DF)
}

# All arguments specified
setMethod("getDistaz", signature(obj="IrisClient", 
                                 latitude="numeric", longitude="numeric",
                                 staLatitude="numeric", staLongitude="numeric"), 
          function(obj, latitude, longitude, staLatitude, staLongitude) 
            getDistaz.IrisClient(obj, latitude, longitude, staLatitude, staLongitude))





