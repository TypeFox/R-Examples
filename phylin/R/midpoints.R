midpoints <-
function(samples, x=1, y=2, sp.name=row.names(samples), all=FALSE) {
    
    #samples - data.frame with list of coordinates
    #x - column index or name of x values or longitudes
    #y - column index or name of y values or latitutes
    # names assumed in row.names

    if (all) {
        sp.combn <- t(combn(sp.name, 2))
    } else {
        #require(geometry)
        if (is.element('geometry', installed.packages()[,1]) )
        {
            tri <- geometry::delaunayn(samples[,c(x,y)])
            sp.combn <- cbind(sp.name[c(tri[,1], tri[,2], tri[,3])], 
                              sp.name[c(tri[,2], tri[,3], tri[,1])])

            sp.combn <- t(apply(sp.combn, 1, sort)) # order so (2,1) = (1,2)
            sp.combn <- sp.combn[!duplicated(sp.combn),] # remove duplicates
        } else {
            stop("Package 'geometry' is needed for Delaunay triangulation.")
        }
    }

    #source sample and target sample
    colnames(sp.combn) <- c('ss', 'ts') 
    nc <- nrow(sp.combn)

    midpoints <- data.frame(sp.combn, x=rep(NA, nc), y=rep(NA, nc))
    m1 <- match(midpoints$ss, sp.name)
    m2 <- match(midpoints$ts, sp.name)
    midpoints$x <- (samples[m1,x] + samples[m2,x]) / 2
    midpoints$y <- (samples[m1,y] + samples[m2,y]) / 2

    midpoints

}
