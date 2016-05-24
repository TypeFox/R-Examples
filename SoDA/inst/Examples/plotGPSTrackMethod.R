plotGPSTrack <- function(x, y,
                         colors = c("blue", "green", "brown"), 
                         fraction = .75, head = .5, ...) {
    elevationColors <- .arrowColors(x@elevation, colors)
    plotGPSArrows(x, fraction, head, col = elevationColors)
}
