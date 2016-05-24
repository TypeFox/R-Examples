# Verified 1.3.18
# Version 5.0
get.Grid.size <-
function(shape, origin.grid, x.res=0.05, y.res=0.05, plot=FALSE) {
        g = get.shape.range(shape)
        p3 = g[3]
        p1 = g[1]
        p2 = g[2]
        p4 = g[4]
        if ( ! missing(origin.grid) ) {
                p3 = origin.grid$lats[origin.grid$lats <= g[3]]
                if ( length(p3) == 0 ) { stop("origin.grid does NOT contain the current shape latitudinal wise. Expand grid.") }
                p3 = p3[length(p3)]
                p1 = origin.grid$longs[origin.grid$longs <= g[1]]
                if ( length(p1) == 0 ) { stop("origin.grid does NOT contain the current shape longitudinal wise. Expand grid.") }
                p1 = p1[length(p1)]
                p2 = origin.grid$longs[origin.grid$longs >= g[2]]
                if ( length(p2) == 0 ) { stop("origin.grid does NOT contain the current shape longitudinal wise. Expand grid.") }
                p2 = p2[1]
                p4 = origin.grid$lats[origin.grid$lats >= g[4]]
                if ( length(p4) == 0 ) { stop("origin.grid does NOT contain the current shape latitudinal wise. Expand grid.") }
                p4 = p4[1]
                if ( (origin.grid$x.res / x.res)%%1 != 0 ) { stop("origin.grid's x.res in NOT a multiple given x.res.") }
                if ( (origin.grid$y.res / y.res)%%1 != 0 ) { stop("origin.grid's y.res in NOT a multiple given y.res.") }
        }
        lats = seq(p3, p4, y.res)
        longs = seq(p1, p2, x.res)
        if ( range(lats)[2] < p4 ) { lats = c(lats, lats[length(lats)]+y.res) }
        if ( range(longs)[2] < p2 ) { longs = c(longs, longs[length(longs)]+x.res) }
        hrow = length(lats)
        hcol = length(longs)
        if ( plot == T ) {
                if ( ! missing(origin.grid) ) { 
                        plot(rep(origin.grid$longs, each=origin.grid$nrow), rep(origin.grid$lats, origin.grid$ncol), asp=1, pch=3, xlab=NA, ylab=NA)
                        points(rep(longs, each=hrow), rep(lats,hcol), pch=4, col='red', cex=0.5)
                } else {
                        plot(rep(longs, each=hrow), rep(lats,hcol), pch=4, col='red', xlab=NA, ylab=NA, asp=1)
                }
                plot(shape, add=TRUE)
                abline(v=g[, 1:2], lty=2)
                abline(h=g[, 3:4], lty=2)
                cat("Grid size: ", hrow, "x", hcol, "\n")
        }
        return(list(ncol=hcol, nrow=hrow, longs=longs, lats=lats, x.res=x.res, y.res=y.res, grid=cbind(rep(longs, each=hrow), rep(lats, hcol))))
}
