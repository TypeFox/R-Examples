# Verified 1.3.18
get.shape.range <-
function(shape) {
        SHP.range = matrix(ncol=4, nrow=length(shape))
        for ( i in 1:length(shape) ) {
                d = slot(shape, "polygons")[[i]]
                SHP.sub = matrix(ncol=4, nrow=length(slot(d, "Polygons")))
                for ( j in 1 : length(slot(d, "Polygons")) ) {
                        d.sub = slot(d, "Polygons")[[j]]
                        d.sub = slot(d.sub, "coords")
                        SHP.sub[j, 1:2] = range(d.sub[,1])
                        SHP.sub[j, 3:4] = range(d.sub[,2])
                }
                d = matrix(apply(SHP.sub, 2, range), ncol=4)
                SHP.range[i, 1:2] = diag(d[1:2,1:2])
                SHP.range[i, 3:4] = diag(d[1:2,3:4])
        }
        d = matrix(apply(SHP.range, 2, range), ncol=2)
        d = matrix(apply(d, 2, range), ncol=4)
        colnames(d) <- c("Long.start", "Long.end", "Lat.start", "Lat.end")
        return(d)
}
