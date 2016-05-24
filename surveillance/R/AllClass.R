# -------------  class sts  ----------------------------------------

.sts <- setClass(
    "sts",
    slots = c(
        epoch = "numeric",  # this slot was called "week" in surveillance < 1.3
        freq = "numeric",
        start = "numeric",
        observed = "matrix",
        state = "matrix",
        alarm = "matrix",
        upperbound = "matrix",
        neighbourhood = "matrix",
        populationFrac = "matrix",
        map = "SpatialPolygons",
        control = "list",
        ## New slots added in version 1.1-2 to handle proportion time series:
        epochAsDate = "logical",
        multinomialTS = "logical"
    ),
    prototype = list(
        start = c(2000, 1), freq = 52,  # historical defaults
        epochAsDate = FALSE, multinomialTS = FALSE
    ),
    validity = function (object) {
        dimObserved <- dim(object@observed)
        namesObserved <- colnames(object@observed)
        ## CAVE: NULL colnames are uncommon but possible,
        ##       e.g., after aggregate(object, by="unit")
        errors <- c(
            if (!isScalar(object@freq) || object@freq <= 0)
                "'freq' must be a single positive number",
            if (length(object@start) != 2)
                "'start' must be of length two: (year, week/month/idx)",
            if (!is.numeric(object@observed))
                "'observed' must be a numeric matrix",
            ## check consistency of slot dimensions wrt dim(observed):
            if (length(object@epoch) != dimObserved[1L])
                "'epoch' must be of length 'nrow(observed)'",
            if (!identical(dim(object@state), dimObserved))
                "'state' must have the same dimensions as 'observed'",
            if (!identical(dim(object@alarm), dimObserved))
                "'alarm' must have the same dimensions as 'observed'",
            if (!identical(dim(object@upperbound), dimObserved))
                "'upperbound' must have the same dimensions as 'observed'",
            if (!identical(dim(object@neighbourhood), dimObserved[c(2L,2L)]))
                "'neighbourhood' must be a square matrix of size 'ncol(observed)'",
            if (!identical(dim(object@populationFrac), dimObserved))
                "'populationFrac' must have the same dimensions as 'observed'",
            ## if a map is provided, it must cover all colnames(observed):
            if (length(object@map) > 0 && # i.e., not the empty prototype
                !all(namesObserved %in% row.names(object@map)))
                "'map' is incomplete; ensure that all(colnames(observed) %in% row.names(map))",
            ## FIXME: invalidate objects with a map but NULL colnames?
            ## check booleans
            if (length(object@epochAsDate) != 1 || is.na(object@epochAsDate))
                "'epochAsDate' must be either TRUE or FALSE",
            if (length(object@multinomialTS) != 1 || is.na(object@multinomialTS))
                "'multinomialTS' must be either TRUE or FALSE"
        )
        if (length(errors) > 0) errors else TRUE
    }
)


######################################################################
# Definition of the stsBP class for backprojections.
######################################################################

setClass("stsBP",
         slots = list(
             ci = "array",
             lambda = "array"
         ),
         contains = "sts")


######################################################################
# Definition of the stsNC class for nowcasts.
######################################################################

setClass("stsNC",
         slots = list(
             reportingTriangle = "matrix",
             predPMF = "list",
             pi = "array",
             truth = "sts",
             delayCDF = "list",
             SR = "array"
         ),
         contains = "sts")
