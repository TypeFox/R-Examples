potency <-
function (object) 
{
    if ((class(object) == "pla") | (class(object) == "plaCRD") | 
        (class(object) == "plaLSDD") | (class(object) == "plaRBD")) {
        fit <- pla.fit(object@data, sampleLabels = object@sampleLabels, 
                       design = object@design, dfAdj = object@dfAdjustment, 
                       dr = object@dilutionRatio, main = object@assayTitle, 
                       alpha = object@alpha, factor = object@factor, show = FALSE, 
                       returnPotencyEstimates = TRUE)
        fit@pheur$KP
    }
    else {
        if (class(object) == "plaFit") {
            object@pheur$KP
        }
    }
}
