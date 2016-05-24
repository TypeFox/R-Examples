# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## sample control
setMethod("length", "VirtualSampleControl", function(x) getK(x))

## sample setup
setMethod("length", "SampleSetup", function(x) length(getIndices(x)))

## contamination control
setMethod("length", "VirtualContControl", function(x) length(getEpsilon(x)))

## NA control
setMethod("length", "VirtualNAControl", function(x) getLength(getNArate(x)))

getLength <- function(x) {
    if(is(x, "numeric")) length(x) 
    else if(is(x, "matrix")) nrow(x)
    else NA  # other classes
}

