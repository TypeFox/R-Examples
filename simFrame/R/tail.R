# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## data control
setMethod("tail", "VirtualDataControl", function(x) x)

## sample control
setMethod("tail", "VirtualSampleControl", function(x) x)

## sample setup
setMethod("tail", "SampleSetup", 
    function(x, k = 6, n = 6, ...) {
        if(!is.numeric(k) || length(k) == 0) k <- 6
        else k <- k[1]
        indices <- tail(getIndices(x), n=k, ...)  # last list components
        indices <- lapply(indices, tail, n=n, ...)  # last elements of components
        setIndices(x, indices)
        call <- match.call(call=sys.call(-1))  # jump back one environment
        setCall(x, call)
        x
    })

## contamination control
setMethod("tail", "VirtualContControl", function(x) x)

## NA control
setMethod("tail", "VirtualNAControl", function(x) x)

## strata information
setMethod("tail", "Strata", 
    function(x, ...) {
        values <- getValues(x)
        indices <- tail(1:length(values), ...)
        values <- values[indices]
        nr <- getNr(x)
        n <- length(indices)
        if(n == 0) split <- replicate(length(nr), integer())
        else {
            split <- split(indices, factor(values, levels=nr))
            split <- unname(split)
        }
        size <- sapply(split, length)
        call <- match.call()
        Strata(values=values, split=split, design=getDesign(x), 
            nr=nr, legend=getLegend(x), size=size, call=call)
    })

## simulation control
setMethod("tail", "SimControl", function(x) x)


# class "SimResults"
setMethod("tail", "SimResults", 
    function(x, ...) {
        values <- tail(getValues(x), ...)
        setValues(x, values)
        call <- match.call()
        setCall(x, call)
        x
    })
