LowUpFromKnots <- function (knots,verbose=FALSE) { ## assumes knots are in samplingSpace, as LowUpfn does
    localLow <- apply(knots, 2, min)
    localUp <- apply(knots, 2, max)
    return(LowUpfn(localLow, localUp, boundstype = "numerical",verbose=verbose))
}
