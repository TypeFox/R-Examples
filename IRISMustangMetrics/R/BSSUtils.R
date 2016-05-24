##
##    Utility functions for communicating with the Backend Storage Service (BSS)
##    of IRIS DMC project MUSTANG.
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
# Method to save a list of Metric objects to disk.
#
# Internet access methods come from the RCurl package which is loaded by the
# IRISSeismic package which is required by the siesmicMetrics package.
################################################################################

saveMetricList <- function(metricList,
                           id=Sys.getpid(), rdata=FALSE) {
  
  name <- paste0("metricList_",id)      
  
  if (rdata) {
    
    # Save the data as .RData file
    filename <- paste(name,"RData",sep=".")
    result <- try( saveReturn <- save(metricList,file=filename),
                   silent=TRUE)
    if (class(result)[1] == "try-error" ) {
      stop(geterrmessage())      
    }
    
  } else {
    
    # Save the data as .xml file
    filename <- paste(name,"xml",sep=".")
    
    # Choose the appropriate data-to-XML conversion
    if (length(metricList) == 1 && class(metricList[[1]]) == "character") { # For latency which comes preformatted as XML
      xml <- metricList[[1]]
    } else if (class(metricList[[1]]) == "SingleValueMetric") {
      xml <- metricList2Xml(metricList)    
    } else if (class(metricList[[1]]) == "SpectrumMetric") {
      xml <- spectrumMetric2Xml(metricList)
    } else if (length(metricList) == 1 && class(metricList[[1]]) == "MultipleTimeValueMetric") {
      xml <- timesMetric2Xml(metricList[[1]])
    } else {
      stop(paste("saveMetricList: metric of class '",class(metricList[[1]]),"' is not recognized.",sep=""))
    }
    
    result <- try( saveReturn <- write(file=filename,xml),
                   silent=TRUE)
    if (class(result)[1] == "try-error" ) {
      stop(geterrmessage())      
    }
    
  }
  
  return(invisible(filename))
  
}

################################################################################
# Methods to parse and simplifiy incoming error messages
################################################################################

convertBssErrors <- function(err_msg) {
  
  # Look for lines like this:
  #
  # <message>...</message>
  
  # TODO:  Use XML parsing to get BSS error messages
  # NOTE:  It seems that the long java dumps are cut off when returned by geterrmessage() so we don't have complete XML
  # NOTE:  That's why we're using simple string matching right now.
  lines <- unlist(stringr::str_split(err_msg,"\n"))
  bssErrorMessages <- lines[stringr::str_detect(lines,"<message>")]
  err_msg <- bssErrorMessages
  err_msg <- stringr::str_trim(stringr::str_replace(stringr::str_replace(err_msg,"</message>",""),"<message>",""))
  return(err_msg)
}


################################################################################
# getBssMetricList method returns a list of SingleValueMetric objects
#
# This method extends the IrisClient object from package "IRISSeismic" so that it 
# can access BSS webservices that are internal to project MUSTANG.
#
# See documentation on Duncan Temple Lang's PDF package for details on the 
# XML parsing functions used in this method.
################################################################################

if (!isGeneric("getBssMetricList")) {
  setGeneric("getBssMetricList", function(obj, network, station, location, channel,
                                          starttime, endtime, metricName, url) {
    standardGeneric("getBssMetricList")
  })
}

getBssMetricList.IrisClient <- function(obj, network, station, location, channel,
                                        starttime, endtime, metricName, url) {
  
  # TODO:  Talk with Andrew C about harmonizing BSS with webservice wildcarding syntax
  # Convert webservice style wildcarding with "*" to BSS equivalent "".
  network <- stringr::str_replace(network, "\\*", "")
  station <- stringr::str_replace(station, "\\*", "")
  location <- stringr::str_replace(location, "\\*", "")
  channel <- stringr::str_replace(channel, "\\*", "")
  
  # Obtain the XML and parse it
  result <- try( measurementsXml <- getMetricsXml(obj,network,station,location,channel,starttime,endtime,metricName,url),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getBssMetricList.IrisClient: URL Not Found:",url))
    } else {
      stop(paste("getBssMetricList.IrisClient:",err_msg))
    } 
    
  }  
  
  # No errors so proceed
  
  doc <- XML::xmlRoot(XML::xmlTreeParse(measurementsXml))
  
  if (length(XML::xmlChildren(doc)) == 0) {
    
    # Empty result -- no measurements found
    metricList <- list()
    
  } else {
    
    # TODO:  Add support for multiple error messages
    # Check for errors
    bssErrors <- XML::xmlElementsByTagName(doc,"error")
    if (length(bssErrors) > 0) {
      msg1 <- XML::xmlElementsByTagName(bssErrors[[1]])[[1]]
      err_msg <- paste("getBssMetricList:",XML::xmlValue(msg1))
      stop(err_msg)
      ###MCRQuit(err_msg,"BSS_QUERY")
    }
    
    if (metricName == "up_down_times") {
      
      # An example returning a single MultipleTimeValueMetric
      #
      # http://service.iris.edu//mustang/measurements/1/query?metric=up_down_times&target=IU.XMAS.10.BHE.B&startafter=2012-07-00T00:00:00@endbefore=2012-07-00T00:00:00
      # <measurements>
      #  <up_down_times net="IU" sta="XMAS" loc="10" chan="BHE" qual="M" target="IU.XMAS.10.BHE.B" start="2012-07-01T00:00:00.000" end="2012-07-02T00:00:00.000">
      #   <t value="2012-07-01T00:00:00.000"/>
      #   <t value="2012-07-01T22:48:36.000"/>
      #   <t value="2012-07-01T22:58:43.000"/>
      #   <t value="2012-07-01T23:23:31.000"/>
      #   <t value="2012-07-01T23:33:39.000"/>
      #   <t value="2012-07-01T23:59:59.000"/>
      #  </up_down_times>
      # </measurements>
      
      # TODO:  support BSS response with more than one MultipleTimeValueMetric
      if (length(XML::xmlChildren(doc)) > 1) {
        err_msg <- paste("getBssMetricList:  more than one MultipleTimeValueMetric in BSS response")
        stop(err_msg)
        ###MCRQuit(err_msg,"BSS_QUERY")        
      }
      
      # Extract measurements attributes as a character matrix with one column for each measurement
      # and rows named: "value"  "net"    "sta"    "loc"    "chan"   "qual"   "target" "start"  "end"
      # NOTE:  The following attributeMatrix should only have a single row
      attributeMatrix <- XML::xmlSApply(doc, XML::xmlAttrs)
      
      # Pre-allocate a list of the proper length
      metricList <- vector("list", ncol(attributeMatrix))
      
      # Convert information into MultipleTimeValueMetrics
      # NOTE:  Should oly be one metric
      for (i in seq(ncol(attributeMatrix))) {
        snclq <- attributeMatrix["target",i]
        starttime <- as.POSIXct(attributeMatrix["start",i], format="%Y-%m-%dT%H:%M:%OS", tz="GMT")
        endtime <- as.POSIXct(attributeMatrix["end",i], format="%Y-%m-%dT%H:%M:%OS", tz="GMT")
        metricName <- colnames(attributeMatrix)[i]
        
        # NOTE:  <t value="..."/> elements only have a single attribute so we don't have to acces it with ["value"]
        metricNodes <- XML::xmlChildren(doc)
        timeStrings <- as.character(XML::xmlSApply(metricNodes[[1]],XML::xmlAttrs))
        upDownTimes <- as.POSIXct(timeStrings, format="%Y-%m-%dT%H:%M:%OS", tz="GMT") 
        
        metricList[[i]] <- new("MultipleTimeValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="up_down_times", values=upDownTimes)        
      }
      
      
    } else {
      
      # An example returning multiple SingleValueMetrics
      #
      # http://service.iris.edu/mustang/measurements/1/query?net=IU&sta=ANMO&loc=00&chan=BH?&timewindow=2012-07-12T00:00:00,2012-07-13T00:00:00&metric=percent_availability
      # <measurements>
      #  <percent_availability value="100.00" net="IU" sta="ANMO" loc="00" chan="BH1" qual="M" target="IU.ANMO.00.BH1.B" start="2012-07-12T00:00:00.000" end="2012-07-13T00:00:00.000"/>
      #  <percent_availability value="100.00" net="IU" sta="ANMO" loc="00" chan="BH2" qual="M" target="IU.ANMO.00.BH2.B" start="2012-07-12T00:00:00.000" end="2012-07-13T00:00:00.000"/>
      #  <percent_availability value="100.00" net="IU" sta="ANMO" loc="00" chan="BHZ" qual="M" target="IU.ANMO.00.BHZ.B" start="2012-07-12T00:00:00.000" end="2012-07-13T00:00:00.000"/>
      #  ...
      # </measurements>
      
      # Extract measurements attributes as a character matrix with one column for each measurement
      # and rows named: "value"  "net"    "sta"    "loc"    "chan"   "qual"   "target" "start"  "end"
      attributeMatrix <- XML::xmlSApply(doc, XML::xmlAttrs)
      
      # Pre-allocate a list of the proper length
      metricList <- vector("list", ncol(attributeMatrix))
      
      # Convert information into SingleValueMetrics
      for (i in seq(ncol(attributeMatrix))) {
        snclq <- attributeMatrix["target",i]
        starttime <- as.POSIXct(attributeMatrix["start",i], format="%Y-%m-%dT%H:%M:%OS", tz="GMT")
        endtime <- as.POSIXct(attributeMatrix["end",i], format="%Y-%m-%dT%H:%M:%OS", tz="GMT")
        metricName <- colnames(attributeMatrix)[i]
        value <- as.numeric(attributeMatrix["value",i])
        metricList[[i]] <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName=metricName, value=value)
      }
      
    }
    
  }
  
  return(metricList)
}

# All arguments specified
setMethod("getBssMetricList", signature(obj="IrisClient", 
                                        network="character", location="character", station="character", 
                                        channel="character",starttime="POSIXct", endtime="POSIXct", 
                                        metricName="character", url="character"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, url) 
            getBssMetricList.IrisClient(obj, network, station, location, channel, starttime, endtime, 
                                        metricName, url))
# url="missing"
setMethod("getBssMetricList", signature(obj="IrisClient", 
                                        network="character", location="character", station="character", 
                                        channel="character",starttime="POSIXct", endtime="POSIXct", 
                                        metricName="character", url="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, url) 
            getBssMetricList.IrisClient(obj, network, station, location, channel,
                                        starttime, endtime, metricName,
                                        url="http://service.iris.edu/mustang/measurements/1/query?"))


################################################################################
# getMetricsXml method returns measurementsXml from the BSS metrics webservice:
#
# http://service.iris.edu/mustang/measurements
#
################################################################################

if (!isGeneric("getMetricsXml")) {
  setGeneric("getMetricsXml", function(obj, network, station, location, channel,
                                       starttime, endtime, metricName, url) {
    standardGeneric("getMetricsXml")
  })
}

getMetricsXml.IrisClient <- function(obj, network, station, location, channel,
                                     starttime, endtime, metricName, url) {
  
  # Assemble URL
  if (network != "") {
    url <- paste(url,"net=",network,sep="")    
  }
  if (station != "") {
    url <- paste(url,"&sta=",station,sep="")
  }
  if (location != "") {
    # NOTE:  Unlike other web services, the measurements service expects "" rather than "--" for the BLANK location designator.
    location <- stringr::str_replace(location,"--","")
    url <- paste(url,"&loc=",location,sep="")
  }
  if (channel != "") {
    url <- paste(url,"&cha=",channel,sep="")
  }
  # TODO:  What sort of starttime/endtime checking should we have?
  url <- paste(url,"&timewindow=",format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),",",
               format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="") # webservice requires "T" format
  
  url <- paste(url,"&metric=",metricName,sep="")
  
  if (obj@debug) {
    write(paste("URL =",url), stdout())
  }
  
  # Get data from BSS measurements webservice
  # NOTE:  RCurl::getURLContent returns a binary objected based on the "resulting HTTP  header's Content-Type field."
  # NOTE:  Use RCurl::getURL to return data as character.
  result <- try( measurementsXml <- RCurl::getURL(url, useragent=obj@useragent),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getMetricsXml.IrisClient: URL Not Found:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getMetricsXml.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getMetricsXml.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  return(measurementsXml)
}

# All arguments specified
setMethod("getMetricsXml", signature(obj="IrisClient", 
                                     network="character", location="character", station="character", 
                                     channel="character",starttime="POSIXct", endtime="POSIXct", 
                                     metricName="character", url="character"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, url) 
            getMetricsXml.IrisClient(obj, network, station, location, channel,
                                     starttime, endtime, metricName, url))
# url="missing"
setMethod("getMetricsXml", signature(obj="IrisClient", 
                                     network="character", location="character", station="character", 
                                     channel="character",starttime="POSIXct", endtime="POSIXct", 
                                     metricName="character", url="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, url) 
            getMetricsXml.IrisClient(obj, network, station, location, channel,
                                     starttime, endtime, metricName,
                                     url="http://service.iris.edu/mustang/measurements/1/query?"))



################################################################################
# getLatencyValuesXml method returns XML from the LatencyValues servlet:
#
# http://www.iris.edu/mustang/latency/LatencyValues?net=IU&stn=ANMO&loc=*&chn=BHZ
#
################################################################################

if (!isGeneric("getLatencyValuesXml")) {
  setGeneric("getLatencyValuesXml", function(obj, network, station, location, channel, url) {
    standardGeneric("getLatencyValuesXml")
  })
}

getLatencyValuesXml.IrisClient <- function(obj, network, station, location, channel, url) {
  
  # Assemble URL
  if (network != "") {
    url <- paste(url,"net=",network,sep="")    
  }
  if (station != "") {
    url <- paste(url,"&sta=",station,sep="")
  }
  if (location != "") {
    # NOTE:  Unlike other web services, the measurements service expects "" rather than "--" for the BLANK location designator.
    location <- stringr::str_replace(location,"--","")
    url <- paste(url,"&loc=",location,sep="")
  }
  if (channel != "") {
    url <- paste(url,"&cha=",channel,sep="")
  }
  
  if (obj@debug) {
    write(paste("URL =",url), stdout())
  }
  
  # Get data from LatencyValues servlet
  # NOTE:  RCurl::getURLContent returns a binary objected based on the "resulting HTTP  header's Content-Type field."
  # NOTE:  Use RCurl::getURL to return data as character.
  result <- try( latencyValuesXml <- RCurl::getURL(url, useragent=obj@useragent),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getLatencyValuesXml.IrisClient: URL Not Found:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getLatencyValuesXml.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getLatencyValuesXml.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  return(latencyValuesXml)
}

# All arguments specified
setMethod("getLatencyValuesXml", signature(obj="IrisClient", 
                                           network="character", location="character", station="character", 
                                           channel="character", url="character"), 
          function(obj, network, station, location, channel, url) 
            getLatencyValuesXml.IrisClient(obj, network, station, location, channel, url))
# url="missing"
setMethod("getLatencyValuesXml", signature(obj="IrisClient", 
                                           network="character", location="character", station="character", 
                                           channel="character", url="missing"), 
          function(obj, network, station, location, channel, url) 
            getLatencyValuesXml.IrisClient(obj, network, station, location, channel,
                                           url="http://www.iris.edu/mustang/latency/LatencyValues?"))


################################################################################
# Simpler access to measurements from MUSTANG
################################################################################

# createBssUrl -----------------------------------------------------------------
#
# Create the URL needed to access data from the MUSTANG BSS
#
# Exaple of getting single valued measurements output from the ''measurements' service:
#
# http://service.iris.edu/mustang/measurements/1/query?net=IU&sta=ANMO&loc=10&cha=BHZ&timewindow=2013-06-01T00:00:00,
#       2013-06-02T00:00:00%20&output=text&metric=sample_mean,sample_min

if (!isGeneric("createBssUrl")) {
  setGeneric("createBssUrl", function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) {
    standardGeneric("createBssUrl")
  })
}

createBssUrl.IrisClient <- function(obj, network, station, location, channel,
                                    starttime, endtime, metricName, constraint, url) {
  
  # Assemble URL
  if (network != "") {
    url <- paste(url,"net=",network,sep="")    
  }
  if (station != "") {
    url <- paste(url,"&sta=",station,sep="")
  }  
  if (location != "") {
    # NOTE:  Unlike other web services, the measurements service expects "" rather than "--" for the BLANK location designator.
    location <- stringr::str_replace(location,"--","")
    url <- paste(url,"&loc=",location,sep="")
  }
  if (channel != "") {
    url <- paste(url,"&cha=",channel,sep="")
  }
  # REC -- allow a new constraint term called 'intersects' which means we don't use the timewindow term
  # intersection means that metric start is less than specified end of time window and
  # metric end is greater than start of time window.
  if (constraint == "intersects") {
    # &start_lt=2012-08-08T00:00:00&end_gt=2012-08-02T00:00:00  
    url <- paste(url,"&start_lt=",
                 format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),"&end_gt=",
                 format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="") # webservice requires "T" format
  } else {
    url <- paste(url,"&timewindow=",
                 format(starttime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),",",
                 format(endtime,"%Y-%m-%dT%H:%M:%OS", tz="GMT"),sep="") # webservice requires "T" format
  }
  
  url <- paste(url,"&output=text",sep="")
  # metricName can be a comma-separated list of metric names
  url <- paste(url,"&metric=",metricName,sep="")
  
  # Add any "constraint" argument verbatim
  if (constraint != "" && constraint != "intersects") {
    url <- paste(url,"&",constraint,sep="")    
  }
  
  if (obj@debug) {
    write(paste("URL =",url), stdout())
  }
  
  return(url)  
  
}

# All arguments specified
setMethod("createBssUrl", signature(obj="IrisClient", 
                                    network="character", location="character", station="character", 
                                    channel="character",starttime="POSIXct", endtime="POSIXct", 
                                    metricName="character", constraint="character", url="character"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) 
            createBssUrl.IrisClient(obj, network, station, location, channel, starttime, endtime, 
                                    metricName, constraint, url))
# url="missing"
setMethod("createBssUrl", signature(obj="IrisClient", 
                                    network="character", location="character", station="character", 
                                    channel="character",starttime="POSIXct", endtime="POSIXct", 
                                    metricName="character", constraint="character", url="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) 
            createBssUrl.IrisClient(obj, network, station, location, channel,
                                    starttime, endtime, metricName, constraint,
                                    "http://service.iris.edu/mustang/measurements/1/query?"))

# constraint="missing"
setMethod("createBssUrl", signature(obj="IrisClient", 
                                    network="character", location="character", station="character", 
                                    channel="character",starttime="POSIXct", endtime="POSIXct", 
                                    metricName="character", constraint="missing", url="character"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) 
            createBssUrl.IrisClient(obj, network, station, location, channel,
                                    starttime, endtime, metricName, "", url))


# constraint="missing", url="missing"
setMethod("createBssUrl", signature(obj="IrisClient", 
                                    network="character", location="character", station="character", 
                                    channel="character",starttime="POSIXct", endtime="POSIXct", 
                                    metricName="character", constraint="missing", url="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) 
            createBssUrl.IrisClient(obj, network, station, location, channel,
                                    starttime, endtime, metricName, "",
                                    "http://service.iris.edu/mustang/measurements/1/query?"))



# getSingleValueMetrics ---------------------------------------------------
#
# NOTE:  metricName can be a comma-separated list of metric names.
#
# Exaple of getting single valued measurements output from the ''measurements' service:
#
# http://service.iris.edu/mustang/measurements/1/query?net=IU&sta=ANMO&loc=10&cha=BHZ&timewindow=2013-06-01T00:00:00,
#       2013-06-02T00:00:00%20&output=text&metric=sample_mean,sample_min
#
# "Sample Mean Metric"
# "value","target","start","end","lddate"
# "-7021.92","IU.ANMO.10.BHZ.B","2013/06/01 00:00:00","2013/06/02 00:00:00","1970/01/01 00:00:00"
# "-7021.92","IU.ANMO.10.BHZ.M","2013/06/01 00:00:00","2013/06/02 00:00:00","2013/09/11 04:28:08.394736"
# 
# "Sample Min Metric"
# "value","target","start","end","lddate"
# "-19622.0","IU.ANMO.10.BHZ.B","2013/06/01 00:00:00","2013/06/02 00:00:00","1970/01/01 00:00:00"
# "-19622.0","IU.ANMO.10.BHZ.M","2013/06/01 00:00:00","2013/06/02 00:00:00","2013/09/11 04:28:08.394736"

if (!isGeneric("getSingleValueMetrics")) {
  setGeneric("getSingleValueMetrics", function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) {
    standardGeneric("getSingleValueMetrics")
  })
}

getSingleValueMetrics.IrisClient <- function(obj, network, station, location, channel,
                                             starttime, endtime, metricName, constraint, url) {
  
  # Create the BSS URL
  url <- createBssUrl(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url)  
  
  # Get data from the measurements webservice
  # NOTE:  RCurl::getURLContent returns a binary objected based on the "resulting HTTP  header's Content-Type field."
  # NOTE:  Use RCurl::getURL to return data as character.
  result <- try( measurementsText <- RCurl::getURL(url, useragent=obj@useragent),
                 silent=TRUE)
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getSingleValueMetrics.IrisClient: URL Not Found:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getSingleValueMetrics.IrisClient: couldn't connect to host"))
    } 
    
  } else {
    
    # Handle error messages coming directly from the BSS
    # NOTE:  Encountered a situation where measurementsText had more than one element
    lines <- unlist(stringr::str_split(measurementsText,"\\n"))
    if ( stringr::str_detect(lines[1],"Exception") ) {
      err_msg <- lines[2]
      stop(paste("getSingleValueMetrics.IrisClient:",err_msg))
    }
    
  }
  
  # No errors so proceed
  
  # Dataframes will be returned in a list
  dataframeList <- list()
  
  # NOTE:  The metrics will not necessarily be returned in the order requested.
  
  # Metric names returned with output=text do not match the requested metric names
  # NOTE:  Copy names from http://service.iris.edu/mustang/measurements/1
  # NOTE:  and do 1 minutes of vim to generate this list.
  convertName <- list("Amplifier Saturation Metric"="amplifier_saturation",
                      "Calibration Signal Metric"="calibration_signal",
                      "Clock Locked Metric"="clock_locked",
                      "Data Latency Metric"="data_latency",
                      "DC Offset Times Metric"="dc_offset_times",
                      "Digital Filter Charging Metric"="digital_filter_charging",
                      "Digitizer Clipping Metric"="digitizer_clipping",
                      "Event Begin Metric"="event_begin",
                      "Event End Metric"="event_end",
                      "Event In Progress Metric"="event_in_progress",
                      "Feed Latency Metric"="feed_latency",
                      "Glitches Metric"="glitches",
                      "Max Gap Metric"="max_gap",
                      "Max LTA/STA Metric"="max_ltasta",
                      "Max Overlap Metric"="max_overlap",
                      "Max STA/LTA Metric"="max_stalta",
                      "Missing Padded Data Metric"="missing_padded_data",
                      "Num Gaps Metric"="num_gaps",
                      "Num Overlaps Metric"="num_overlaps",
                      "Num Spikes Metric"="num_spikes",
                      "Percent Availability Metric"="percent_availability",
                      "PSD Metric"="psd",
                      "Sample Max Metric"="sample_max",
                      "Sample Mean Metric"="sample_mean",
                      "Sample Median Metric"="sample_median",
                      "Sample Min Metric"="sample_min",
                      "Sample RMS"="sample_rms",
                      "Sample SNR"="sample_snr",
                      "Spikes Metric"="spikes",
                      "Station Completeness Metric"="station_completeness",
                      "Station Up Down Times Metric"="station_up_down_times",
                      "Suspect Time Tag Metric"="suspect_time_tag",
                      "Telemetry Sync Error Metric"="telemetry_sync_error",
                      "Timing Correction Metric"="timing_correction",
                      "Timing Quality Metric"="timing_quality",
                      "Total Latency Metric"="total_latency",
                      "Up Down Times Metric"="up_down_times")
  
  # Break the text into chunks separated by "\n\n".
  # NOTE:  stringr::str_split uses extended regular expressions and '\' needs to be escaped
  chunks <- unlist(stringr::str_split(measurementsText,"\\n\\n"))
  
  # length of chunks represents the number of metrics represented
  for (i in seq(length(chunks))) {
    
    # Create a dataframe from the text
    DF <- utils::read.csv(skip=1, header=TRUE, stringsAsFactors=FALSE, text=chunks[i])
    
    # Get metric name from first line of chunk
    lines <- unlist(stringr::str_split(chunks[i],"\\n"))
    bssName <- stringr::str_replace_all(lines[1],'"','')
    measurementName <- convertName[[bssName]]
    
    # REC
    # if there is no measurementName found in the convertName list, let's take the metric
    # name that is the header line for this chunk and convert it to the BSS metric name.
    if (length(measurementName) == 0) {
      # split text components on space
      splitNameList <- strsplit(bssName,"[[:blank:]]")
      # convert text to lower case
      toLowerNameList <- unlist(lapply(splitNameList,tolower))
      # remove the last term ('metric')
      shorterName <- toLowerNameList[-(length(toLowerNameList))]
      # join text with underscores
      measurementName <- paste(shorterName,collapse="_")
    }
    
    # Convert from BSS DF names to 'seismic' package standard naming
    names <- names(DF)
    for (i in seq(length(names))) {
      if (names[i] == 'value') {
        names[i] <- measurementName
      } else if (names[i] == 'target') {
        names[i] <- "snclq"
      } else if (names[i] == 'start') {
        names[i] <- "starttime"
      } else if (names[i] == 'end') {
        names[i] <- "endtime"
      } else if (names[i] == 'lddate') {
        names[i] <- "loadtime"
      }
    }
    names(DF) <- names
    
    # Convert time strings
    DF$starttime <- as.POSIXct(DF$starttime, "%Y/%m/%d %H:%M:%OS", tz="GMT")
    DF$endtime <- as.POSIXct(DF$endtime, "%Y/%m/%d %H:%M:%OS", tz="GMT")    
    DF$loadtime <- as.POSIXct(DF$loadtime, "%Y/%m/%d %H:%M:%OS", tz="GMT")
    
    # NOTE:  The database was originally populated with a version of this package
    # NOTE:  that always assigned quality to be 'B'. Later versions obtained the
    # NOTE:  quality from the miniSEED packet (typically 'M').  Because of this
    # NOTE:  it is possible to have duplicate entries that only differ in the Q
    # NOTE:  part of their snclq.  To avoid double counting, we need to remove 
    # NOTE:  duplicates when the only difference is 'B'/'M'.
    
    # Reorder rows to be in temporal order by reverse loadtime -- most recent first
    DF <- DF[order(DF$loadtime,decreasing=TRUE),]
    
    # Create a uniqueId that does not use quality by first removing the Q part of N.S.L.C.Q
    sncl <- stringr::str_replace(DF$snclq,"\\.[A-Z]$","")
    uniqueId <- paste(sncl,format(DF$starttime,"%Y%m%d"),format(DF$endtime,"%Y%m%d"),sep='.')
    
    # Remove any rows which share a duplicate uniqueId (older loadtime, i.e. quality='B')
    DF <- DF[!duplicated(uniqueId),]
    
    # Reorder rows to be in temporal order by starttime
    DF <- DF[order(DF$starttime),]
    measurementName_DF <- paste(measurementName,"DF",sep="_")
    dataframeList[[measurementName_DF]] <- DF
  }
  
  # NOTE:  Previously, we just returned the list of dataframes
  ###return(dataframeList)  
  
  # Now convert dataframeList into a 'tidy' dataframe appropriate for use with ggplot2
  
  # This function will be applied to each dataframe in DFList
  # Adds a "metricName" column and "value" column
  addMetricName <- function(df) {
    metric <- names(df)[[1]]
    df$metricName <- metric
    n <- names(df)
    n[[1]] <- "value"
    names(df) <- n
    return(df)
  }
  
  # Apply function to our list and bind the rows together to make one big dataframe
  # NOTE:  http://www.r-bloggers.com/the-rbinding-race-for-vs-do-call-vs-rbind-fill/
  TidyDF <- do.call("rbind", lapply(dataframeList, addMetricName))
  
  # Move last column to first 
  TidyDF <- TidyDF[,c(6,1,2,3,4,5)]
  
  # Remove ugly rownames
  row.names(TidyDF) <- NULL
  
  return(TidyDF)
  
}

# All arguments specified
setMethod("getSingleValueMetrics", signature(obj="IrisClient", 
                                             network="character", location="character", station="character", 
                                             channel="character",starttime="POSIXct", endtime="POSIXct", 
                                             metricName="character", constraint="character", url="character"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) 
            getSingleValueMetrics.IrisClient(obj, network, station, location, channel, starttime, endtime, 
                                             metricName, constraint, url))
# url="missing"
setMethod("getSingleValueMetrics", signature(obj="IrisClient", 
                                             network="character", location="character", station="character", 
                                             channel="character",starttime="POSIXct", endtime="POSIXct", 
                                             metricName="character", constraint="character", url="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) 
            getSingleValueMetrics.IrisClient(obj, network, station, location, channel,
                                             starttime, endtime, metricName, constraint,
                                             "http://service.iris.edu/mustang/measurements/1/query?"))

# constraint="missing"
setMethod("getSingleValueMetrics", signature(obj="IrisClient", 
                                             network="character", location="character", station="character", 
                                             channel="character",starttime="POSIXct", endtime="POSIXct", 
                                             metricName="character", constraint="missing", url="character"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) 
            getSingleValueMetrics.IrisClient(obj, network, station, location, channel,
                                             starttime, endtime, metricName, "", url))


# constraint="missing", url="missing"
setMethod("getSingleValueMetrics", signature(obj="IrisClient", 
                                             network="character", location="character", station="character", 
                                             channel="character",starttime="POSIXct", endtime="POSIXct", 
                                             metricName="character", constraint="missing", url="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, metricName, constraint, url) 
            getSingleValueMetrics.IrisClient(obj, network, station, location, channel,
                                             starttime, endtime, metricName, "",
                                             "http://service.iris.edu/mustang/measurements/1/query?"))

# getPsdMetrics -----------------------------------------------------------
#
# Exaple of getting PSD output from the 'measurements' service
#
# http://service.iris.edu/mustang/measurements/1/query?metric=psd&net=IU&sta=ANMO&cha=BHZ&loc=00&output=tex&timewindow=2013-05-01T00:00:00,2013-05-01T01:00:00
#
# "PSD Metric"
# "target","start","end","f","a","p"
# "IU.ANMO.00.BHZ.M","2013/05/01 00:30:00","2013/05/01 01:30:00","10","-15.254799999999999","0"
# "IU.ANMO.00.BHZ.M","2013/05/01 00:30:00","2013/05/01 01:30:00","9.1700400000000002","-14.224399999999999","0"
# "IU.ANMO.00.BHZ.M","2013/05/01 00:30:00","2013/05/01 01:30:00","8.4089600000000004","-13.209199999999999","0"
# ... ... through all frequencies
# ... repeated for each hour long chunk

if (!isGeneric("getPsdMetrics")) {
  setGeneric("getPsdMetrics", function(obj, network, station, location, channel, starttime, endtime, url) {
    standardGeneric("getPsdMetrics")
  })
}

getPsdMetrics.IrisClient <- function(obj, network, station, location, channel, starttime, endtime, url) {
  
  # TODO:  Sanity check for arguments
  
  # Create the BSS URL
  url <- createBssUrl(obj, network, station, location, channel, starttime, endtime, metricName="psd", url=url)  
  
  # Read the data directly from the URL with utils::read.csv
  
  colNames <- c("target","starttime","endtime","frequency","amplitude","phase")
  
  # REC -- modify this call to make use of RCurl
  #result <- try( DF <- utils::read.csv(url, header=FALSE, skip=2, col.names=colNames, stringsAsFactors=FALSE),
  #               silent=TRUE )
  result <- try( DF <- utils::read.csv(textConnection(RCurl::getURL(url,useragent=obj@useragent)), header=FALSE, skip=2, col.names=colNames, stringsAsFactors=FALSE),
                 silent=TRUE )
  
  # Handle error response
  if (class(result) == "try-error" ) {
    
    err_msg <- geterrmessage()
    if (stringr::str_detect(err_msg,"Not Found")) {
      stop(paste("getPsdMetrics.IrisClient: URL Not Found:",url))
    } else if (stringr::str_detect(err_msg,"couldn't connect to host")) {
      stop(paste("getPsdMetrics.IrisClient: couldn't connect to host"))
    } else {
      stop(paste("getPsdMetrics.IrisClient:",err_msg))
    } 
    
  }
  
  # No errors so proceed
  
  # NOTE:  Using 'output=text' returns non-ISO datetimes instead of the ISO 
  # NOTE:  datetimes returned when 'output=xml'.
  
  # Convert time strings
  DF$starttime <- as.POSIXct(DF$starttime, "%Y/%m/%d %H:%M:%OS", tz="GMT")
  DF$endtime <- as.POSIXct(DF$endtime, "%Y/%m/%d %H:%M:%OS", tz="GMT")
  
  # Convert phase from integer to numeric
  DF$phase <- as.numeric(DF$phase)
  
  return(DF)
}

# All arguments specified
setMethod("getPsdMetrics", signature(obj="IrisClient", 
                                     network="character", location="character", station="character", 
                                     channel="character", starttime="POSIXct", endtime="POSIXct", 
                                     url="character"), 
          function(obj, network, station, location, channel, starttime, endtime, url) 
            getPsdMetrics.IrisClient(obj, network, station, location, channel, starttime, endtime, url))
# url="missing"
setMethod("getPsdMetrics", signature(obj="IrisClient", 
                                     network="character", location="character", station="character", 
                                     channel="character", starttime="POSIXct", endtime="POSIXct", 
                                     url="missing"), 
          function(obj, network, station, location, channel, starttime, endtime, url) 
            getPsdMetrics.IrisClient(obj, network, station, location, channel, starttime, endtime,
                                     url="http://service.iris.edu/mustang/measurements/1/query?"))


