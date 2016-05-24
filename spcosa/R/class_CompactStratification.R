setClass(
    Class = "CompactStratification",
    representation = representation(
        cells = "SpatialPixels",
        stratumId = "integer",
        centroids = "SpatialPoints",
        mssd = "numeric"
    ),
    contains = "Stratification"
)
