# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## data control
setMethod("head", "VirtualDataControl", function(x) x)

## sample control
setMethod("head", "VirtualSampleControl", function(x) x)

## sample setup
setMethod("head", "SampleSetup", 
    function(x, k = 6, n = 6, ...) {
        if(!is.numeric(k) || length(k) == 0) k <- 6
        else k <- k[1]
        indices <- head(getIndices(x), n=k, ...)  # first list components
        indices <- lapply(indices, head, n=n, ...)  # first elements of components
        setIndices(x, indices)
        call <- match.call(call=sys.call(-1))  # jump back one environment
        setCall(x, call)
        x
    })

## contamination control
setMethod("head", "VirtualContControl", function(x) x)

## NA control
setMethod("head", "VirtualNAControl", function(x) x)

## strata information
setMethod("head", "Strata", 
    function(x, ...) {
        values <- head(getValues(x), ...)
        nr <- getNr(x)
        n <- length(values)
        if(n == 0) split <- replicate(length(nr), integer())
        else {
            split <- split(1:length(values), factor(values, levels=nr))
            split <- unname(split)
        }
        size <- sapply(split, length)
        call <- match.call()
        Strata(values=values, split=split, design=getDesign(x), 
            nr=nr, legend=getLegend(x), size=size, call=call)
    })

## simulation control
setMethod("head", "SimControl", function(x) x)

# class "SimResults"
setMethod("head", "SimResults", 
    function(x, ...) {
        values <- head(getValues(x), ...)
        setValues(x, values)
        call <- match.call()
        setCall(x, call)
        x
    })
