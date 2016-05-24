## This function is taken and slightly modified from the example of ?identify.
## It uses identify to select points, and recolor the selected points.

pickCenter <- function(x, y = NULL, n = length(x),
                       pch = 19, col = "red", cex = 1.2,
                       labels = seq_along(x), labelcex = 1, ...)
    {
        xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
        sel <- rep(FALSE, length(x)); res <- integer(0)
        N <- 1
        while(sum(sel) < n) {
            ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
            if(!length(ans)) break
            ans <- which(!sel)[ans]
            points(x[ans], y[ans], pch = pch,
                   cex = cex, col = col[(N - 1) %% length(col) + 1])
            text(x[ans], y[ans], labels = ans, pos = 1, cex = labelcex)
            sel[ans] <- TRUE
            res <- c(res, ans)
            N <- N + 1
        }
        if(sum(sel) > length(col))
            warning(paste0("Number of clusters (",sum(sel),
                           ") > number of colors(",length(col),
                           "). Colors recycled."))
        invisible(res)
    }
