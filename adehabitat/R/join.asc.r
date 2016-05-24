"join.asc" <- function(pts, x)
{
    ## Verifications
    if (!inherits(x, "asc")) stop("non convenient data")

    ## coordinates of the limits of the pixels
    xy<-getXYcoords(x)
    xy$x<-xy$x+attr(x, "cellsize")/2
    xy$x<-c(xy$x, xy$x[1]-attr(x, "cellsize")/2)
    xy$y<-xy$y+attr(x, "cellsize")/2
    xy$y<-c(xy$y, xy$y[1]-attr(x, "cellsize")/2)

    ## cuts the points and stores into two vectors of
    ## indices (one for rows and one for columns)
    xf<-as.numeric(cut(pts[,1], xy$x))
    yf<-as.numeric(cut(pts[,2], xy$y))

    ## In case of a factor map, stores the levels
    fact<-0
    if (attr(x, "type")=="factor")
      ct<-attr(x, "levels")

    ## For each point
    for (i in 1:nrow(pts)) {

        ## identifies the value of the map
        if (attr(x, "type")=="numeric") {
            u<-x[xf[i],yf[i]]
            fact[i]<-u
        }
        if (attr(x, "type")=="factor") {
            u<-x[xf[i],yf[i]]
            tmp<-ct[u]
            if (length(tmp)==1) {
                fact[i]<-tmp
            } else {
                fact[i]<-NA
            }
        }
    }
    ## output
    if (attr(x, "type")=="factor") fact<-factor(fact)
    return(fact)
}

