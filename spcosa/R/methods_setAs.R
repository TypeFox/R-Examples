setAs(
    from = "CompactStratification",
    to = "data.frame",
    def = function(from) {
        as(as(from, "SpatialPixelsDataFrame"), "data.frame")
    }
)



setAs(
    from = "CompactStratification",
    to = "SpatialPixels",
    def = function(from) {
        tmp <- from@cells
        gridded(tmp) <- FALSE
        suppressWarnings(as(tmp, "SpatialPixels"))
    }
)



setAs(
    from = "CompactStratification",
    to = "SpatialPixelsDataFrame",
    def = function(from) {
        suppressWarnings(
            SpatialPixelsDataFrame(
                points = as(from, "SpatialPixels"), # type cast is needed for class "SpatialPixelsDataFrame"
                data = data.frame(stratumId = from@stratumId)
            )
        )
    }
)



setAs(
    from = "SamplingPattern",
    to = "data.frame",
    def = function(from) {
        as(from@sample, "data.frame")
    }
)



setAs(
    from = "SamplingPatternRandomComposite",
    to = "data.frame",
    def = function(from) {
        as(as(from, "SpatialPointsDataFrame"), "data.frame")
    }
)



setAs(
    from = "SamplingPatternRandomComposite",
    to = "SpatialPointsDataFrame",
    def = function(from) {
        SpatialPointsDataFrame(
            coords = as(from, "SpatialPoints"),
            data = data.frame(composite = from@composite)
        )
    }
)



setAs(
    from = "SamplingPattern",
    to = "SpatialPoints",
    def = function(from) {
        from@sample
    }
)
