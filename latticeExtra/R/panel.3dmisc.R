

## a panel function for cloud that draws "3d bar charts" which
## shouldn't be used except in unusual circumstances


panel.3dbars <-
    function(x, y, z,
             rot.mat = diag(4), distance,
             xbase = 1, ybase = 1,
             xlim, xlim.scaled,
             ylim, ylim.scaled,
             zlim, zlim.scaled,
             zero.scaled,
             col = "black",
             lty = 1, lwd = 1,
             alpha = 1,
             ...,
             col.facet = "white",
             alpha.facet = 1)
{
    n <- length(z)
    col <- rep(col, length = n)
    col.facet <- rep(col.facet, length = n)
    alpha <- rep(alpha, length = n)
    alpha.facet <- rep(alpha.facet, length = n)
    lty <- rep(lty, length = n)
    lwd <- rep(lwd, length = n)
    id <-
        ((x >= xlim.scaled[1]) & (x <= xlim.scaled[2]) &
         (y >= ylim.scaled[1]) & (y <= ylim.scaled[2]) &
         !is.na(x) & !is.na(y) & !is.na(z))
    m <- ltransform3dto3d(rbind(x, y, 0), rot.mat, distance)
    ord <- sort.list(m[3,])
    ord <- ord[id[ord]]
    zero.scaled <-
        if (zero.scaled < zlim.scaled[1]) zlim.scaled[1]
        else if (zero.scaled > zlim.scaled[2]) zlim.scaled[2]
        else zero.scaled
    inRange <- function(x, lim)
    {
        rng <- range(x, finite = TRUE)
        rng[1] >= min(lim) && rng[2] <= max(lim)
    }
    ## draw bars one by one
    for (i in ord)
    {
        ## print(i)
        xbase.scaled <- diff(xlim.scaled) * xbase / diff(xlim)
        ybase.scaled <- diff(ylim.scaled) * ybase / diff(ylim)
        zz.sides <- matrix(c(zero.scaled, z[i]), 2, 1)[, rep(1, 5)]
        xx.sides <-
            c(x[i], x[i]) + xbase.scaled * 0.5 *
                rbind(c(-1, 1, 1, -1, -1), c(-1, 1, 1, -1, -1))
        yy.sides <-
            c(y[i], y[i]) + ybase.scaled * 0.5 *
                rbind(c(-1, -1, 1, 1, -1), c(-1, -1, 1, 1, -1))
        zz.top <- matrix(z[i], 2, 2)
        xx.top <- 
            c(x[i], x[i]) + xbase.scaled * 0.5 *
                rbind(c(-1, 1), c(-1, 1))
        yy.top <-
            c(y[i], y[i]) + ybase.scaled * 0.5 *
                rbind(c(-1, -1), c(1, 1))

        zz <- cbind(zz.sides, c(NA, NA), zz.top)
        xx <- cbind(xx.sides, c(NA, NA), xx.top)
        yy <- cbind(yy.sides, c(NA, NA), yy.top)
        ## str(list(xx, yy, zz))
        if (inRange(xx, xlim.scaled) &&
            inRange(yy, ylim.scaled) &&
            inRange(zz, zlim.scaled))
        {
            panel.3dwire(xx, yy, zz,
                         rot.mat = rot.mat, distance = distance,
                         xlim = xlim, xlim.scaled = xlim.scaled,
                         ylim = ylim, ylim.scaled = ylim.scaled,
                         zlim = zlim, zlim.scaled = zlim.scaled,
                         col = col[i], lty = lty[i], lwd = lwd[i],
                         alpha = alpha[i],
                         ...,
                         at = c(0, 1), # dummy
                         col.regions = col.facet[i],
                         alpha.regions = alpha.facet[i])
        }
    }
}




## panel.3dpolygon <-
##     function(x, y, z, rot.mat = diag(4), distance,
##              type = 'p',
##              xlim.scaled,
##              ylim.scaled,
##              zlim.scaled,
##              zero.scaled,
##              col = "white", 
##              border = "black", 
##              lty = 1, lwd = 1,
##              min.sides = 3,
##              ...,
##              subscripts = TRUE)
## {
##     m <- ltransform3dto3d(rbind(x, y, z), rot.mat, distance)
##     ## ord <- sort.list(m[3,])
##     n <- ncol(m)
##     w <- which(is.na(x) | is.na(y))
##     id.lengths <- diff(c(0, w, n))
##     cum.lengths <- c(0, cumsum(id.lengths))

##     idlist <-
##         lapply(seq_along(id.lengths),
##                function(i) {
##                    ind <- seq_len(id.lengths[i]) + cum.lengths[i]
##                    ind[-id.lengths[i]]
##                })

##     ord <-
##         order(sapply(idlist,
##                      function(ind) {
##                          min(m[3, ind])
##                      }))

##     for (ind in idlist[ord])
##     {
##         if (length(ind) >= min.sides)
##             panel.polygon(x = m[1, ind], y = m[2, ind],
##                           col = col, border = border)
##     }

## }




panel.3dpolygon <-
    function(x, y, z, rot.mat = diag(4), distance,
             xlim.scaled,
             ylim.scaled,
             zlim.scaled,
             zero.scaled,
             col = "white", 
             border = "black", 
             ## min.sides = 3,
             font, fontface, ## gpar() doesn't like these
             ...)
{
    if (all(is.na(x) | is.na(y) | is.na(z))) return()
    border <- 
        if (all(is.na(border)))
            "transparent"
        else if (is.logical(border))
        {
            if (border) "black"
            else "transparent"
        }
        else border
    m <- ltransform3dto3d(rbind(x, y, z), rot.mat, distance)
    ## ord <- sort.list(m[3,])
    n <- ncol(m)
    w <- which(is.na(x) | is.na(y))
    id.lengths <- diff(c(0, w, n))

    ## need to reorder multiple polygons by some measure of "average" depth

    id.long <- rep(seq_along(id.lengths), id.lengths)
    ord.depth <- order(tapply(m[3,], id.long, min, na.rm = TRUE))
    id.ordered <- ord.depth[id.long]

    ord.long <- order(id.ordered)
    
    grid.polygon(x = m[1, ord.long], y = m[2, ord.long],
                 id = id.ordered[ord.long],
                 default.units = "native",
                 gp =
                 gpar(fill = col,
                      col = border,
                      ...))

##     print(data.frame(x = m[1, ord.long],
##                      y = m[2, ord.long],
##                      id = id.ordered[ord.long]))

    return()
}




panel.3dtext <-
    function(x, y, z, labels = seq_along(x),
             rot.mat = diag(4), distance, ...)
{
    if (all(is.na(x) | is.na(y) | is.na(z))) return()
    m <- ltransform3dto3d(rbind(x, y, z), rot.mat, distance)
    ord <- sort.list(m[3,])
    panel.text(x = m[1, ord], y = m[2, ord], labels = labels, ...)
}



## d <- data.frame(x = rnorm(10),
##                 y = rnorm(10),
##                 z = rnorm(10))
## rownames(d) <- letters[1:10]

## cloud(z ~ x * y, d, panel.3d.cloud = panel.3dtext)
## cloud(z ~ x * y, d, panel.3d.cloud = panel.3dtext,
##       labels = rownames(d), col = "red")

## ## for multipanel plots

## cloud(z ~ x * y, d,
##       panel.3d.cloud = function(..., subscripts) {
##           panel.3dtext(..., labels = rownames(d)[subscripts])
##       })

