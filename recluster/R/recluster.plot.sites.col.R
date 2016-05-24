recluster.plot.sites.col<-function (long, lat, mat, cext = 0.3, cex = 1, cex.axis = 0.7, 
    cex.lab = 0.8, text = FALSE, pch=21, add = FALSE,...) 
{
    if(!add) {plot(long, lat, col = rgb(mat[, 3], mat[, 4], mat[, 
        5], maxColorValue = 255), cex = cex, cex.axis = cex.axis, 
        cex.lab = cex.lab, pch = pch,...)}
    if(add) {points(long, lat, col = rgb(mat[, 3], mat[, 4], mat[, 
        5], maxColorValue = 255), cex = cex, cex.axis = cex.axis, 
        cex.lab = cex.lab, pch = pch,...)}
}
