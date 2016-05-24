# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## data control
setMethod("summary", "VirtualDataControl", function(object) object)

## sample control
setMethod("summary", "VirtualSampleControl", function(object) object)

## sample setup
setMethod("summary", "SampleSetup", 
    function(object) {
        size <- sapply(getIndices(object), length, USE.NAMES=FALSE)
        if(length(size) == 0) size <- numeric()  # in case of empty list
        SummarySampleSetup(size=size)
    })

## contamination control
setMethod("summary", "VirtualContControl", function(object) object)

## NA control
setMethod("summary", "VirtualNAControl", function(object) object)

## strata information
setMethod("summary", "Strata", 
    function(object) data.frame(getLegend(object), Size=getSize(object)))

## simulation control
setMethod("summary", "SimControl", function(object) object)

## simulation results
# summary of data.frame is of class "table", 
# thus no extra class "summary.SimResults" necessary
#setMethod("summary", "SimResults", 
#    function(object, ...) summary(getValues(object), ...))
setMethod("summary", "SimResults", 
    function(object, ...) {
        values <- getValues(object)
#        values$Run <- NULL
        values$Run <- as.factor(values$Run)
        if(!is.null(values$Sample)) values$Sample <- as.factor(values$Sample)
        if(!is.null(values$Rep)) values$Rep <- as.factor(values$Rep)
        if(!is.null(values$Epsilon)) values$Epsilon <- as.factor(values$Epsilon)
        if(!is.null(values$NArate)) values$NArate <- as.factor(values$NArate)
        summary(values, ...)
    })
