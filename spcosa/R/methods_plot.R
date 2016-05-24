setMethod(
    f = "plot",
    signature = signature(
        x = "CompactStratification",
        y = "missing"
    ),
    definition = function(x, y, ...) {
        
        # convert to "SpatialPixelsDataFrame"
        suppressWarnings( # suppress warning 'grid has empty column/rows in dimension 2'
            x <- as(x, "SpatialPixelsDataFrame")
        )
        
        # extract grid topology
        gridTopology <- getGridTopology(x)
        cellSize <- getCellSize(x)
        
        # simple row and columwise edge detection algorithm
        M <- as(x, "matrix")
        i <- 2:nrow(M)
        j <- 2:ncol(M)
        Bc <- rbind(   (M[i,, drop = FALSE] - M[i-1,, drop = FALSE]), 0)
        Br <- cbind(0, (M[,j, drop = FALSE] - M[,j-1, drop = FALSE]))
        Bc[Bc != 0] <- 1
        Br[Br != 0] <- 1
        Bc <- SpatialGridDataFrame(grid = gridTopology, data = data.frame(edge = as.integer(Bc)))
        Br <- SpatialGridDataFrame(grid = gridTopology, data = data.frame(edge = as.integer(Br)))
        Bc <- as(Bc, "data.frame")
        Br <- as(Br, "data.frame")
        Bc <- Bc[Bc$edge == 1, coordnames(x)]
        Br <- Br[Br$edge == 1, coordnames(x)]
        tmp <- Br
        Br$x <- tmp[[1]] - 0.5 * cellSize[1]
        Br$y <- tmp[[2]] + 0.5 * cellSize[2]
        Br$xend <- Br$x + cellSize[1]
        Br$yend <- Br$y
        tmp <- Bc
        Bc$x <- tmp[[1]] + 0.5 * cellSize[1]
        Bc$y <- tmp[[2]] - 0.5 * cellSize[2]
        Bc$xend <- Bc$x
        Bc$yend <- Bc$y + cellSize[2]
        
        # plot map and strata boundaries
        s <- coordinates(as(x, "SpatialPixels"))
        sNames <- colnames(s)
        if (is(sNames, "NULL")) {
            sNames <- c("", "")
        }
        colnames(s) <- c("s1", "s2")
        p <- ggplot() +
            geom_raster(
                data = as.data.frame(s),
                mapping = aes_string(x = "s1", y = "s2"),
                fill   = rgb(0.5, 0.8, 0.5, 1.0)
            ) +
            coord_fixed() +
            scale_x_continuous(name = sNames[1]) +
            scale_y_continuous(name = sNames[2])
        if (nrow(Br) > 0) {
            p <- p +
                geom_segment(
                    data = Br,
                    mapping = aes_string(x = "x", y = "y", xend = "xend", yend = "yend"),
                    colour = rgb(0.8, 1.0, 0.8)
                )
        }
        if (nrow(Bc) > 0) {
            p <- p +
                geom_segment(
                    data = Bc,
                    mapping = aes_string(x = "x", y = "y", xend = "xend", yend = "yend"),
                    colour = rgb(0.8, 1.0, 0.8)
                )
        }
        p + theme(legend.position = "none")
    }
)



setMethod(
    f = "plot",
    signature = signature(
        x = "CompactStratification",
        y = "SamplingPattern"
    ),
    definition = function(x, y, ...) {
        s <- coordinates(as(y, "SpatialPoints"))
        colnames(s) <- c("s1", "s2")
        plot(x) +
            geom_point(
                data = as.data.frame(s),
                mapping = aes_string(
                    x = "s1",
                    y = "s2"
                ),
                colour = "black",
                alpha = 0.4
            )
    }
)



setMethod(
    f = "plot",
    signature = signature(
        x = "CompactStratification",
        y = "SamplingPatternPriorPoints"
    ),
    definition = function(x, y, ...) {
        s <- coordinates(as(y, "SpatialPoints"))
        colnames(s) <- c("s1", "s2")
        plot(x) +
            geom_point(
                data = data.frame(
                    s,
                    type = ifelse(
                        y@isPriorPoint,
                        "prior point",
                        "new point"
                    )
                ),
                mapping = aes_string(
                    x = "s1",
                    y = "s2",
                    shape = "type"
                ),
                colour = "black",
                alpha = 0.4
            ) +
            theme(legend.position = "right")
    }
)



setMethod(
    f = "plot",
    signature = signature(
        x = "CompactStratification",
        y = "SamplingPatternRandomComposite"
    ),
    definition = function(x, y, ...) {
        s <- coordinates(as(y, "SpatialPoints"))
        colnames(s) <- c("s1", "s2")
        plot(x) +
            geom_point(
                data = data.frame(
                    s,
                    composite = factor(y@composite)
                ),
                mapping = aes_string(
                    x = "s1",
                    y = "s2",
                    shape = "composite"
                ),
                colour = "black",
                alpha = 0.4
            ) +
            theme(legend.position = "right")
    }
)



setMethod(
    f = "plot",
    signature = signature(
        x = "SamplingPattern",
        y = "missing"
    ),
    definition = function(x, y, ...) {
        s <- coordinates(as(x, "SpatialPoints"))
        sNames <- colnames(s)
        if (is(sNames, "NULL")) {
            sNames <- c("", "")
        }
        colnames(s) <- c("s1", "s2")
        ggplot() +
            geom_point(
                data = as.data.frame(s),
                mapping = aes_string(x = "s1", y = "s2"),
                colour = "black",
                alpha = 0.4
            ) +
            coord_fixed() +
            scale_x_continuous(name = sNames[1]) +
            scale_y_continuous(name = sNames[2])
    }
)



setMethod(
    f = "plot",
    signature = signature(
        x = "SamplingPatternPriorPoints",
        y = "missing"
    ),
    definition = function(x, y, ...) {
        s <- coordinates(as(x, "SpatialPoints"))
        sNames <- colnames(s)
        if (is(sNames, "NULL")) {
            sNames <- c("", "")
        }
        colnames(s) <- c("s1", "s2")
        ggplot() +
            geom_point(
                data = data.frame(
                    s,
                    type = ifelse(
                        x@isPriorPoint,
                        "prior point",
                        "new point"
                    )
                ),
                mapping = aes_string(
                    x = "s1",
                    y = "s2",
                    shape = "type"
                ),
                colour = "black",
                alpha = 0.4
            ) +
            coord_fixed() +
            scale_x_continuous(name = sNames[1]) +
            scale_y_continuous(name = sNames[2]) +
            theme(legend.position = "right")
    }
)



setMethod(
    f = "plot",
    signature = signature(
        x = "SamplingPatternRandomComposite",
        y = "missing"
    ),
    definition = function(x, y, ...) {
        s <- coordinates(as(x, "SpatialPoints"))
        sNames <- colnames(s)
        if (is(sNames, "NULL")) {
            sNames <- c("", "")
        }
        colnames(s) <- c("s1", "s2")
        ggplot() +
            geom_point(
                data = data.frame(
                    s,
                    composite = factor(x@composite)
                ),
                mapping = aes_string(
                    x = "s1",
                    y = "s2",
                    shape = "composite"
                ),
                colour = "black",
                alpha = 0.4
            ) +
            coord_fixed() +
            scale_x_continuous(name = sNames[1]) +
            scale_y_continuous(name = sNames[2])
    }
)
