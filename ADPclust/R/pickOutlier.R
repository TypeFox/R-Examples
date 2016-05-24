## This function is taken and slightly modified from the example of ?identify.
## Here is the original description:
## A function to use identify to select points, and overplot the
## points with another symbol as they are selected

pickOutlier <- function(x, y = NULL, n = length(x), pch = 19, 
                        col = "black", cex = 1.2, labels = seq_along(x), labelcex = 0.8, ...)
    {
        xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
        sel <- rep(FALSE, length(x)); res <- integer(0)
        N <- 1
        while(sum(sel) < n) {
            ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
            if(!length(ans)) break
            ans <- which(!sel)[ans]
            points(x[ans], y[ans], pch = pch,
                    cex = cex, col = col)

            text(x[ans], y[ans], labels = ans, pos = 3, cex = labelcex)
            sel[ans] <- TRUE
            res <- c(res, ans)
            N <- N + 1
        }
        invisible(res)
    }
