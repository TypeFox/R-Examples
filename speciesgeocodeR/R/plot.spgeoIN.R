plot.spgeoIN <- function(x, ...) {
    xmax <- min(max(x$species_coordinates[, 1]) + 2, 180)
    xmin <- max(min(x$species_coordinates[, 1]) - 2, -180)
    ymax <- min(max(x$species_coordinates[, 2]) + 2, 90)
    ymin <- max(min(x$species_coordinates[, 2]) - 2, -90)
    difx <- sqrt(xmax^2 + xmin^2)
    dify <- sqrt(ymax^2 + ymin^2)
    if (difx > 90) {
        xmax <- min(xmax + 10, 180)
        xmin <- max(xmin - 10, -180)
        ymax <- min(ymax + 10, 90)
        ymin <- max(ymin - 10, -90)
    }
    map("world", xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    axis(1)
    axis(2)
    box("plot")
    plot(x$polygons, col = "grey60", border = "grey40", add = T)
    points(x$species_coordinates[, 1], x$species_coordinates[, 2], cex = 0.7, pch = 3, col = "blue")
    par(mar = c(4,4,4,4))
} 
