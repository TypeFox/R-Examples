##    S4 classes for handling metrics derived from seismic traces.
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
# R classes for a Metric obj.
#
# This class is part of the MUSTANG project at IRIS DMC.
#
# The Metric class stores generated metrics and allows their output
# as XML.
#
#<measurements>
#  <date start='2012-02-10T00:00:00.000' end='2012-02-10T09:20:00.000'>
#    <target snclq='A.B.C1.D.E'>
#      <example value='1.0'/>
#    </target>
#    <target snclq='A.B.C2.D.E'>
#      <example value='2.0'/>
#    </target>
#    <target snclq='A.B.C3.D.E'>
#      <example value='3.0'/>
#    </target>
#  </date>
#  <date start='2012-02-10T09:20:00.000' end='2012-02-10T18:40:00.000'>
#    <target snclq='A.B.C1.D.E'>
#      <example value='1.1'/>
#    </target>
#    <target snclq='A.B.C2.D.E'>
#      <example value='2.1'/>
#    </target>
#    <target snclq='A.B.C3.D.E'>
#      <example value='3.1'/>
#    </target>
#  </date>
#</measurements>


################################################################################
# Class SingleValueMetric
#
# The basic Metric object contains metadata associated with a value.
#
#   snclq       -- station.network.channel.location.quality identifier
#   starttime   -- starttime for the trace used in this metric
#   endtime     -- endtime for the trace used in this metric
#   metricName  -- name of this metric
#   valueName   -- name of the parameter value stored in the BSS
#   value       -- numeric value of this metric
#   valueString -- string representation of the value for this metric
#   quality_flag -- numeric flag identifying issues in the quality of the calculated metric
#   quality_flagString -- string representation of quality_flag
#   attributeName -- name of an optional, additional attribute
#   attributeValueString -- string represetnation of the attribute value
################################################################################

setClass("SingleValueMetric", 
         # typed slots (aka attributes) for class Trace
         representation(snclq = "character",
                        starttime = "POSIXct",
                        endtime = "POSIXct",
                        metricName = "character",
                        valueName = "character",
                        value = "numeric",
                        valueString = "character",
                        quality_flag = "numeric",
                        quality_flagString = "character",
                        attributeName = "character",
                        attributeValueString = "character"),
         # default values for slots
         prototype(snclq = "",
                   starttime = as.POSIXct("1970-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
                   endtime = as.POSIXct("1970-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
                   metricName = "",
                   valueName = "value",
                   valueString = "NULL",
                   quality_flag = -9,
                   quality_flagString = "-9",
                   attributeName = "",
                   attributeValueString = "")
)

# initialze method
setMethod("initialize", "SingleValueMetric",
          function(.Object, snclq, starttime, endtime, metricName, value, quality_flag,
                   attributeName, attributeValueString, ...) {
              
            .Object@snclq <- snclq
            .Object@starttime <- starttime
            .Object@endtime <- endtime
            .Object@metricName <- metricName
            .Object@value <- value
            
            # Set the quality_flag if it is passed in
            if (!missing(quality_flag)) {
              .Object@quality_flag <- quality_flag
            }
            
            # Set the attributeName if it is passed in
            if (!missing(attributeName)) {
              .Object@attributeName <- attributeName
            }
            
            # Set the attributeValueString if it is passed in
            if (!missing(attributeValueString)) {
              .Object@attributeValueString <- attributeValueString
            }
            
            # metrics whose values are integer
            integerMetrics <- c("num_gaps",
                                "num_glitches",
                                "num_overlaps",
                                "num_pings",
                                "num_spikes",
                                "num_outliers",
                                # State Of Health flags
                                "calibration_signal",
                                "timing_correction",
                                "event_begin",
                                "event_end",
                                "event_in_progress",
                                "clock_locked",
                                "amplifier_saturation",
                                "digitizer_clipping",
                                "spikes",
                                "glitches",
                                "missing_padded_data",
                                "telemetry_sync_error",
                                "digital_filter_charging",
                                "suspect_time_tag")
            
            # Set the valueString
            if (metricName == "example") {
              .Object@valueString <- sprintf("%.2f",value)
            } else if (metricName %in% integerMetrics) {
              .Object@valueString <- sprintf("%d",value)
            } else {
              # default formatting 
              .Object@valueString <- sprintf("%.3f",value)
            }
            
            # Convert missing value from R style to BSS style
            .Object@valueString <- stringr::str_replace(.Object@valueString,"NA","NULL")
            
            # Set the quality_flagString using BSS style missing value
            .Object@quality_flagString <- as.character(.Object@quality_flag)
            
            return(.Object)
          }
)

# show method is called by show() and print() ----------------------------------

show.SingleValueMetric <- function(object) {
  cat ("SingleValueMetric \n")
  cat ("  metric:        " , object@metricName, "\n", sep="")
  cat ("  snclq:         " , object@snclq, "\n", sep="")
  cat ("  starttime:     " , format(object@starttime), "\n", sep="")
  cat ("  endtime:       " , format(object@endtime), "\n", sep="")
  cat ("  value:         " , object@valueString, "\n", sep="")
  for (i in seq(length(object@attributeName))) {
    name <- object@attributeName[i]
    valString <- object@attributeValueString[i]
    n <- 14 - stringr::str_length(name)
    n <- ifelse (n > 0, n, 2)
    if (name != "") {
      cat ("  ",name,":", rep(" ",n), valString, "\n", sep="")      
    }
  } 
}
# NOTE:  method signature must match generic signature for 'show' with argument: 'object'
setMethod("show", signature(object="SingleValueMetric"), function(object) show.SingleValueMetric(object))


################################################################################
# Class MultipleValueMetric
#
# The MultipleValueMetric object contains metadata associated with a list of values.
#
#   snclq         -- station.network.channel.location.quality identifier
#   starttime     -- starttime for the trace used in this metric
#   endtime       -- endtime for the trace used in this metric
#   metricName    -- name of this metric
#   elementName   -- name of the XML element containing the parameter value stored in the BSS
#   values        -- vector of values associated with this metric
#   valueStrings -- vector of string representations of the values associated with this metric
#   quality_flag  -- numeric flag identifying issues in the quality of the calculated metric
#   quality_flagString -- string representation of quality_flag
#
################################################################################

setClass("MultipleValueMetric", 
         # typed slots (aka attributes) for class Trace
         representation(snclq = "character",
                        starttime = "POSIXct",
                        endtime = "POSIXct",
                        metricName = "character",
                        elementName = "character",
                        values = "numeric",
                        valueStrings = "character",
                        quality_flag = "numeric",
                        quality_flagString = "character"),
         # default values for slots
         prototype(snclq = "",
                   starttime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
                   endtime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
                   metricName = "",
                   elementName = "x",
                   valueStrings = "c()",
                   quality_flag = -9,
                   quality_flagString = "-9")
)

# initialze method
setMethod("initialize", "MultipleValueMetric",
          function(.Object, snclq, starttime, endtime, metricName, values, quality_flag, ...) {
            
            .Object@snclq <- snclq
            .Object@starttime <- starttime
            .Object@endtime <- endtime
            .Object@metricName <- metricName
            .Object@values <- values
            
            # metrics whose values are integer
            integerMetrics <- c()
                                
            # Set the valueString
            if (metricName %in% integerMetrics) {
              .Object@valueStrings <- sprintf("%d",values)
            } else {
              # default formatting 
              .Object@valueStrings <- sprintf("%.3f",values)
            }
            
            # Convert missing value from R style to BSS style
            .Object@valueStrings <- stringr::str_replace(.Object@valueStrings,"NA","NULL")
            
            # Set the quality_flag if it is passed in
            if (!missing(quality_flag)) {
              .Object@quality_flag <- quality_flag
            }
            
            # Set the quality_flagString using BSS style missing value
            .Object@quality_flagString <- as.character(.Object@quality_flag)
            
            return(.Object)
          }
)

# show method is called by show() and print() ----------------------------------

show.MultipleValueMetric <- function(object) {
  cat ("MultipleValueMetric \n")
  cat ("  metric:        " , object@metricName, "\n")
  cat ("  snclq:         " , object@snclq, "\n")
  cat ("  starttime:     " , format(object@starttime), "\n")
  cat ("  endtime:       " , format(object@endtime), "\n")
  cat ("  values:        " , object@valueStrings, "\n")
}
# NOTE:  method signature must match generic signature for 'show' with argument: 'object'
setMethod("show", signature(object="MultipleValueMetric"), function(object) show.MultipleValueMetric(object))


################################################################################
# Class SpectrumMetric
#
# The SpectrumMetric object contains data associated with a discrete spectrum.
#
#   snclq         -- station.network.channel.location.quality identifier
#   starttime     -- starttime for the trace used in this metric
#   endtime       -- endtime for the trace used in this metric
#   metricName    -- name of this metric
#   elementName   -- name of the XML element containing the parameter value stored in the BSS
#   freqs         -- vector of spectrum frequencies
#   freqStrings   -- vector of string representations of the spectrum frequencies
#   amps          -- vector of spectrum amplitudes
#   ampStrings    -- vector of string representations of the spectrum amplitudes
#   phases        -- vector of spectrum phases
#   phaseStrings  -- vector of string representations of the spectrum phases
#   quality_flag  -- numeric flag identifying issues in the quality of the calculated metric
#   quality_flagString -- string representation of quality_flag
#
################################################################################

setClass("SpectrumMetric", 
         # typed slots (aka attributes) for class Trace
         representation(snclq = "character",
                        starttime = "POSIXct",
                        endtime = "POSIXct",
                        metricName = "character",
                        elementName = "character",
                        freqs = "numeric",
                        freqStrings = "character",
                        amps = "numeric",
                        ampStrings = "character",
                        phases = "numeric",
                        phaseStrings = "character",
                        quality_flag = "numeric",
                        quality_flagString = "character"),
         # default values for slots
         prototype(snclq = "",
                   starttime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
                   endtime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
                   metricName = "",
                   elementName = "spectrum",
                   freqStrings = "c()",
                   ampStrings = "c()",
                   phaseStrings = "c()",
                   quality_flag = -9,
                   quality_flagString = "-9")
)

# initialze method
setMethod("initialize", "SpectrumMetric",
          function(.Object, snclq, starttime, endtime, metricName, freqs, amps, phases, quality_flag, ...) {
            
            .Object@snclq <- snclq
            .Object@starttime <- starttime
            .Object@endtime <- endtime
            .Object@metricName <- metricName
            .Object@freqs <- freqs
            .Object@amps <- amps
            .Object@phases <- phases
            
            # Set the Strings
            .Object@freqStrings <- sprintf("%g",freqs)
            .Object@ampStrings <- sprintf("%g",amps)
            .Object@phaseStrings <- sprintf("%g",phases)
            
            # Convert missing value from R style to BSS style
            .Object@freqStrings <- stringr::str_replace(.Object@freqStrings,"NA","NULL")
            .Object@ampStrings <- stringr::str_replace(.Object@ampStrings,"NA","NULL")
            .Object@phaseStrings <- stringr::str_replace(.Object@phaseStrings,"NA","NULL")
            
            # Set the quality_flag if it is passed in
            if (!missing(quality_flag)) {
              .Object@quality_flag <- quality_flag
            }
            
            # Set the quality_flagString using BSS style missing value
            .Object@quality_flagString <- as.character(.Object@quality_flag)
            
            return(.Object)
          }
)

# show method is called by show() and print() ----------------------------------

show.SpectrumMetric <- function(object) {
  cat ("SpectrumMetric \n")
  cat ("  metric:        " , object@metricName, "\n")
  cat ("  snclq:         " , object@snclq, "\n")
  cat ("  starttime:     " , format(object@starttime), "\n")
  cat ("  endtime:       " , format(object@endtime), "\n")
  cat ("  freqs:         " , object@freqStrings, "\n")
  cat ("  amps:          " , object@ampStrings, "\n")
  cat ("  phases:        " , object@phaseStrings, "\n")
}
# NOTE:  method signature must match generic signature for 'show' with argument: 'object'
setMethod("show", signature(object="SpectrumMetric"), function(object) show.SpectrumMetric(object))

################################################################################
# Class MultipleTimeValueMetric
#
# The MultipleTimeValueMetric object contains metadata associated with a list of POSIXct values.
#
#   snclq         -- station.network.channel.location.quality identifier
#   starttime     -- starttime for the trace used in this metric
#   endtime       -- endtime for the trace used in this metric
#   metricName    -- name of this metric
#   elementName   -- name of the XML element containing the parameter value stored in the BSS
#   values        -- vector of POSIXct values associated with this metric
#   valueStrings -- vector of string representations of the values associated with this metric
#   quality_flag  -- numeric flag identifying issues in the quality of the calculated metric
#   quality_flagString -- string representation of quality_flag
#
# Here is the first example of what a multi-value metric will looklike:
#
# <measurements>
#   <date start='2012-02-10T00:00:00.000' end='2012-02-10T09:20:00.000'>
#     <target snclq='N.S.L.C1.Q'>
#       <up_down_times>
#         <t value="2012-02-10T00:00:00.000"/>
#         <t value="2012-02-10T00:01:00.000"/>
#         <t value="2012-02-10T00:02:00.000"/>
#         <t value="2012-02-10T00:03:00.000"/>
#       </up_down_times>
#     </target>
#   </date>
# </measurements>
#
################################################################################

setClass("MultipleTimeValueMetric", 
         # typed slots (aka attributes) for class Trace
         representation(snclq = "character",
                        starttime = "POSIXct",
                        endtime = "POSIXct",
                        metricName = "character",
                        elementName = "character",
                        values = "POSIXct",
                        valueStrings = "character",
                        quality_flag = "numeric",
                        quality_flagString = "character"),
         # default values for slots
         prototype(snclq = "",
                   starttime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
                   endtime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
                   metricName = "",
                   elementName = "t",
                   valueStrings = "c()",
                   quality_flag = -9,
                   quality_flagString = "-9")
)

# initialze method
setMethod("initialize", "MultipleTimeValueMetric",
          function(.Object, snclq, starttime, endtime, metricName, values, quality_flag, ...) {
            
            .Object@snclq <- snclq
            .Object@starttime <- starttime
            .Object@endtime <- endtime
            .Object@metricName <- metricName
            .Object@values <- values
            
            # Set the quality_flag if it is passed in
            if (!missing(quality_flag)) {
              .Object@quality_flag <- quality_flag
            }
            
            # Convert values to strings
            .Object@valueStrings <- format(values, format="%Y-%m-%dT%H:%M:%0S", tz="GMT")
            
            # Set the quality_flagString using BSS style missing value
            .Object@quality_flagString <- as.character(.Object@quality_flag)
            
            return(.Object)
          }
)

# show method is called by show() and print() ----------------------------------

show.MultipleTimeValueMetric <- function(object) {
  cat ("MultipleTimeValueMetric \n")
  cat ("  metric:        " , object@metricName, "\n")
  cat ("  snclq:         " , object@snclq, "\n")
  cat ("  starttime:     " , format(object@starttime), "\n")
  cat ("  endtime:       " , format(object@endtime), "\n")
  cat ("  values:        " , object@valueStrings, "\n")
}
# NOTE:  method signature must match generic signature for 'show' with argument: 'object'
setMethod("show", signature(object="MultipleTimeValueMetric"), function(object) show.MultipleTimeValueMetric(object))


################################################################################
# Function to convert a list of SingleValueMetric objects into a list of
# dataframes, one per metricName.
################################################################################

metricList2DFList <- function(metricList) {
  
  # Extract attributes from the list of Metrics
  # NOTE:  use 'sapply' to return a vector as opposed to the list returned by 'lapply'
  snclq <- sapply(metricList, slot, 'snclq')
  starttime <- as.numeric(sapply(metricList, slot, 'starttime'))
  endtime <- as.numeric(sapply(metricList, slot, 'endtime'))
  metricName <- sapply(metricList, slot, 'metricName')
  value <- sapply(metricList, slot, 'value')
  quality_flagString <- sapply(metricList, slot, 'quality_flagString')
  
  # NOTE:  attributeName and attributeValueString can be vectors and are therefore
  # NOTE:  returned as a list of character vectors.
  attributeNameList <- lapply(metricList, slot, 'attributeName')
  attributeValueStringList <- lapply(metricList, slot, 'attributeValueString')
  
  # Create a dataframe from the vectors created above
  df <- as.data.frame(cbind(snclq,starttime,endtime,metricName,value,quality_flagString),
                      stringsAsFactors=FALSE)
  
  # Now add attributes 
  # NOTE:  As we are creating a dataframe, every column must be represented in every row.
  # NOTE:  If singleValueMetrics have different attributes, we will need to include all
  # NOTE:  of them as columns and then insert NA's into those rows where they are not found.
  attNames <- unique(unlist(attributeNameList))
  for (name in attNames) {
    if (name != "") {
      df[[name]] <- rep("",length(snclq))
    }
  }
  
  # Now fill in attribute values where thay are defined
  for (i in seq(length(metricList))) {
    names <- metricList[[i]]@attributeName
    valStrings <- metricList[[i]]@attributeValueString
    for (j in seq(length(names))) {
      name <- names[j]
      valString <- valStrings[j]
      if (name != "") {
        df[[name]][i] <- valString
      }
    }
  }
    
   
  # At this point, all columns are of type "character"
  
  
  dfList <- list()
  
  for (metric in unique(df$metricName)) {
    # Pull out data associated with a single metric
    dfName <- paste0(metric,'_DF')
    sub <- subset(df, metricName == metric)
    sub[[metric]] <- sub$value
    sub$quality_flag <- as.numeric(sub$quality_flagString)
    # Convert columns to appropriate type
    sub$starttime <- as.POSIXct(as.numeric(sub$starttime), origin="1970-01-01", tz="GMT")
    sub$endtime <- as.POSIXct(as.numeric(sub$endtime), origin="1970-01-01", tz="GMT")
    # Get rid of unneeded columns
    sub$metricName <- NULL
    sub$value <- NULL
    sub$quality_flagString <- NULL
    # Get rid of attributes with no values
    for (name in attNames) {
      if ( !any(sub[[name]] != "") ) {
        sub[[name]] <- NULL
      }
    }
    # Add this subset dataframe to the list
    dfList[[dfName]] <- sub
  }
  
  return(dfList)
}

  
  
################################################################################
# Function to convert a list of SingleValueMetric objects into the XML expected by
# the DMC data loader.
################################################################################

metricList2Xml <- function(metricList) {

  returnString <- "<measurements>"
  
  # Extract attributes from the list of Metrics
  # NOTE:  use 'sapply' to return a vector as opposed to the list returned by 'lapply'
  snclq <- sapply(metricList, slot, 'snclq')
  starttime <- as.numeric(sapply(metricList, slot, 'starttime'))
  endtime <- as.numeric(sapply(metricList, slot, 'endtime'))
  metricName <- sapply(metricList, slot, 'metricName')
  valueName <- sapply(metricList, slot, 'valueName')
  valueString <- sapply(metricList, slot, 'valueString')
  quality_flagString <- sapply(metricList, slot, 'quality_flagString')
  
  # NOTE:  attributeName and attributeValueString can be vectors and are therefore
  # NOTE:  returned as a list of character vectors.
  attributeNameList <- lapply(metricList, slot, 'attributeName')
  attributeValueStringList <- lapply(metricList, slot, 'attributeValueString')
  
  # Create a dataframe from the vectors created above
  df <- as.data.frame(cbind(snclq,starttime,endtime,metricName,valueName,valueString,quality_flagString),
                      stringsAsFactors=FALSE)
  
  # Convert columns to appropriate type
  df$starttime <- as.POSIXct(as.numeric(df$starttime), origin="1970-01-01", tz="GMT")
  df$endtime <- as.POSIXct(as.numeric(df$endtime), origin="1970-01-01", tz="GMT")
  
  # Now add attributes 
  # NOTE:  As we are creating a dataframe, every column must be represented in every row.
  # NOTE:  If singleValueMetrics have different attributes, we will need to include all
  # NOTE:  of them as columns and then insert NA's into those rows where they are not found.
  attNames <- unique(unlist(attributeNameList))
  for (name in attNames) {
    if (name != "") {
      df[[name]] <- rep("",length(snclq))
    }
  }
  
  # Now fill in attribute values where thay are defined
  for (i in seq(length(metricList))) {
    names <- metricList[[i]]@attributeName
    valStrings <- metricList[[i]]@attributeValueString
    for (j in seq(length(names))) {
      name <- names[j]
      valString <- valStrings[j]
      if (name != "") {
        df[[name]][i] <- valString
      }
    }
  }  
  
  # Add a dateRange string so that we can find all unique date ranges
  df$dateRange <- paste(df$starttime,df$endtime)
  
  # Organize the metrics by dateRange as in the XML example above
  dateRanges <- sort(unique(df$dateRange))  
  
  for (dateRange in dateRanges) {
    
    # Create a subset of the dataframe that only includes the current dateRange
    dateSubset <- df[df$dateRange == dateRange,]
    
    # Create the string for this dateRange
    # NOTE:  This subset may have multiple rows so we only use dates from the first row
	# REC -- correct date formatting to ISO
    #startString <- as.character(dateSubset$starttime[1])
	startString <- format(dateSubset$starttime[1],"%Y-%m-%dT%H:%M:%OS")
    #endString <- as.character(dateSubset$endtime[1])
	endString <- format(dateSubset$endtime[1],"%Y-%m-%dT%H:%M:%OS")
    dateString <- paste("<date start='",startString,"' end='",endString,"'>",sep="")
    returnString <- paste(returnString,dateString,sep="")
    
    # Organize this subset by snclq
    snclqs <- sort(unique(dateSubset$snclq))
    for (id in snclqs) {
      
      target <- subset(dateSubset, snclq == id)
      
      # Create the string for this snclq
      targetString <- paste("<target snclq='",id,"'>",sep="")
      returnString <- paste(returnString,targetString,sep="")
      
      for (i in seq(nrow(target))) {
        
        # Create the valueString
        valueString <- paste(target$valueName[i],"='",target$valueString[i],"'",sep="")
        
        # Create attributerStrings as needed
        attributeString <- ''
        for (name in attNames) {
          if (name != "") {
            valString <- target[[name]][i]
            if (valString != "") {
              attributeString <- paste(attributeString," ",name,"='",valString,"'",sep="")
            }            
          }
        }
        
        # Create quality_flagString if needed
        if (target$quality_flagString[i] == -9) {
          qualityString <- ''
        } else {
          qualityString <- paste("quality_flag='",target$quality_flagString[i],"'",sep="")
        }
        
        # Create the metricString
        metricString <- paste("<",target$metricName[i]," ",
                              valueString," ",
                              qualityString," ",
                              attributeString," ",
                              "/>",sep="")

#         if (target$quality_flagString[i] == "-9") {
#           # If the quality_flag was not assigned during creation of the metric, don't include it in the XML
#           metricString <- paste("<",target$metricName[i]," ",
#                                 target$valueName[i],"='",target$valueString[i],"'",
#                                 attributeString,
#                                 "/>", sep="")
#           
#         } else {
#           # If the quality_flag was assigned, add it to the XML
#           metricString <- paste("<",target$metricName[i]," ",
#                                 target$valueName[i],"='",target$valueString[i],"' ",
#                                 "quality_flag='",target$quality_flagString[i],"'",
#                                 attributeString,
#                                 "/>" ,sep="")
#         }
        
        returnString <- paste(returnString,metricString,sep="")
        
      } # end of "target" loop for each metric
      
      returnString <- paste(returnString,"</target>",sep="")      
    } # end of "snclq" loop
    
    returnString <- paste(returnString,"</date>",sep="")
  } # end of "date" loop
  
  returnString <- paste(returnString,"</measurements>",sep="")
  
  return(returnString)
}


################################################################################
# Function to convert a MultilpeTimeValueMetric into the XML expected by the 
# DMC data loader.
#
# Here is the first example of what a multi-value metric will looklike:
#
# <measurements>
#   <date start='2012-02-10T00:00:00.000' end='2012-02-10T09:20:00.000'>
#     <target snclq='N.S.L.C1.Q'>
#       <up_down_times>
#         <t value="2012-02-10T00:00:00.000"/>
#         <t value="2012-02-10T00:01:00.000"/>
#         <t value="2012-02-10T00:02:00.000"/>
#         <t value="2012-02-10T00:03:00.000"/>
#       </up_down_times>
#     </target>
#   </date>
# </measurements>
#
################################################################################

timesMetric2Xml <- function(metric) {
  
  returnString <- "<measurements>"
  
  # Add the <date>
  startString <- format(metric@starttime, format="%Y-%m-%dT%H:%M:%OS")
  endString <- format(metric@endtime, format="%Y-%m-%dT%H:%M:%OS")
  dateString <- paste("<date start='",startString,"' end='",endString,"'>",sep="")
  returnString <- paste(returnString,dateString,sep="")
  
  # Add the <target>
  targetString <- paste("<target snclq='",metric@snclq,"'>",sep="")
  returnString <- paste(returnString,targetString,sep="")
  
  # Add the <~metric_name~>
  metricString <- paste("<",metric@metricName,">",sep="")
  returnString <- paste(returnString,metricString,sep="")
  
  # Add the individual elements if any exist
  if (length(metric@valueStrings) > 0) {
    elementsString <- paste("<",metric@elementName," value='",metric@valueStrings,"'/>",sep="",collapse="")
    returnString <- paste(returnString,elementsString,sep="")
  }
  
  # Close up the tags  
  metricCloseString <- paste("</",metric@metricName,">",sep="")
  returnString <- paste(returnString,metricCloseString,sep="")
  returnString <- paste(returnString,"</target>",sep="")      
  returnString <- paste(returnString,"</date>",sep="")
  returnString <- paste(returnString,"</measurements>",sep="")

  return(returnString)
}


################################################################################
# Function to convert a SpectrumMetric into the XML expected by the 
# DMC data loader.
#
# Here is the first example of what a multi-value metric will looklike:
#
# <measurements>
#   <date start='2012-02-10T00:00:00.000' end='2012-02-10T09:20:00.000'>
#     <target snclq='N.S.L.C1.Q'>
#       <psd>
#         <spectrum f="0.01" a="1.045234e4"/>
#         ...
#       </psd>
#     </target>
#   </date>
# </measurements>
#
################################################################################

spectrumMetric2Xml <- function(metricList) {
  
  returnString <- ""
  
  for (metric in metricList) {

    # Add the <date>
    startString <- format(metric@starttime, format="%Y-%m-%dT%H:%M:%OS")
    endString <- format(metric@endtime, format="%Y-%m-%dT%H:%M:%OS")
    dateString <- paste("<date start='",startString,"' end='",endString,"'>",sep="")
    returnString <- paste(returnString,dateString,sep="")
    
    # Add the <target>
    targetString <- paste("<target snclq='",metric@snclq,"'>",sep="")
    returnString <- paste(returnString,targetString,sep="")
    
    # Add the <~metric_name~>
    metricString <- paste("<",metric@metricName,">",sep="")
    returnString <- paste(returnString,metricString,sep="")
    
    # Add the individual elements
    elementsString <- paste("<", metric@elementName,
                            " f='", metric@freqStrings,
                            "' a='", metric@ampStrings,
                            "' p='", metric@phaseStrings,
                            "'/>", sep="", collapse="")
    returnString <- paste(returnString,elementsString,sep="")
    
    # Close up the tags  
    metricCloseString <- paste("</",metric@metricName,">",sep="")
    returnString <- paste(returnString,metricCloseString,sep="")
    returnString <- paste(returnString,"</target>",sep="")      
    returnString <- paste(returnString,"</date>",sep="")
    
  }
  
  returnString <- paste("<measurements>",returnString,"</measurements>",sep="")
  
  return(returnString)
}
