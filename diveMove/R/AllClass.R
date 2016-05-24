
setClass("TDR",
         representation=representation(file="character", dtime="numeric",
             time="POSIXct", depth="numeric", concurrentData="data.frame"),
         prototype=prototype(concurrentData=data.frame()),
         validity=function(object) {
             if (length(object@time) != length(object@depth)) {
                 return("depth and time must have equal lengths")
             }
             time.diffs <- diff(unclass(object@time))
             if (any(time.diffs < 0)) {
                 return("time stamps must be in increasing order")
             }
             if (any(time.diffs == 0)) {
                 return("time stamps must not contain duplicate values")
             }
             ccDataN <- nrow(object@concurrentData)
             if (ccDataN > 0 && ccDataN != length(object@time)) {
                 mes <- paste("concurrentData must have the same number of rows",
                              "as there are time stamps")
                 return(mes)
             }
             if (!slot(object, "dtime")) return("dtime cannot be missing")
             return(TRUE)
         })

.speedNames <- c("velocity", "speed")
setClass("TDRspeed", contains="TDR",
         validity=function(object) {
             ccData <- object@concurrentData
             ccDataNames <- names(ccData)
             speedCol <- ccDataNames %in% .speedNames
             if (length(ccDataNames[speedCol]) != 1) {
                 return("speed is not available in concurrentData slot")
             } else if (!is.numeric(ccData[, speedCol])) {
                 return("speed must be of class numeric")
             }
             return(TRUE)
         })

setClass("TDRcalibrate",
         representation=representation(call="call", tdr="TDR",
           gross.activity="list", dive.activity="data.frame",
           dive.phases="factor", dive.models="list", dry.thr="numeric",
           wet.thr="numeric", dive.thr="numeric", speed.calib.coefs="numeric"),
         prototype=prototype(speed.calib.coefs=c(0, 1)),
         validity=function(object) {
             ndives <- max(object@dive.activity$dive.id, na.rm=TRUE)
             dml <- slot(object, "dive.models")
             if (length(slot(object, "dry.thr")) > 1) {
                 return("dry.thr must be a single number")
             }
             if (length(slot(object, "wet.thr")) > 1) {
                 return("wet.thr must be a single number")
             }
             if (length(slot(object, "dive.thr")) > 1) {
                 return("dive.thr must be a single number")
             }
             if (length(slot(object, "speed.calib.coefs")) != 2) {
                 return("speed.calib.coefs must be a length-2 vector")
             }
             if (length(dml) != ndives) {
                 return("All dives must have a corresponding dive model")
             }
             if (! all(sapply(dml, is, "diveModel"))) {
                 return("All elements of dive.models must be class diveModel")
             }
             return(TRUE)
         })

setOldClass("smooth.spline")
setClass("diveModel",
         representation=representation(label.matrix="matrix",
           dive.spline="smooth.spline", spline.deriv="list",
           descent.crit="numeric", ascent.crit="numeric",
           descent.crit.rate="numeric", ascent.crit.rate="numeric"),
         validity=function(object) {
             if (length(slot(object, "descent.crit")) > 1) {
                 return("descent.crit must be a single number")
             }
             if (length(slot(object, "ascent.crit")) > 1) {
                 return("ascent.crit must be a single number")
             }
             if (length(slot(object, "descent.crit.rate")) > 1) {
                 return("descent.crit.rate must be a single number")
             }
             if (length(slot(object, "ascent.crit.rate")) > 1) {
                 return("ascent.crit.rate must be a single number")
             }
             return(TRUE)
         })

setOldClass("nls")                      # For bout methods
