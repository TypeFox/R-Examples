#### Generating functions for subclasses of interpolRisk
OMSRRisk <- function(samplesize = 100) new("OMSRRisk",type=".OMSE", samplesize = samplesize)
MBRRisk <- function(samplesize = 100)  new("MBRRisk",type=".MBRE", samplesize = samplesize)
RMXRRisk <- function(samplesize = 100)  new("RMXRRisk",type=".RMXE", samplesize = samplesize)
setMethod("samplesize","interpolRisk", function(object)object@samplesize)
setReplaceMethod("samplesize","interpolRisk", function(object, value){
       object@samplesize <- value; object})