.arcp <- function(xy)
{
    if (nrow(xy) < 3)
        return(0);
    x.segmat <- cbind(xy, rbind(xy[2:nrow(xy), ],
                                xy[1, ]))
    abs(sum(x.segmat[,1] * x.segmat[,4] - x.segmat[,3]
            * x.segmat[,2])) / 2
}



.arcpspdf <- function(spdf)
{
    lar <- unlist(lapply(polygons(spdf)@polygons,
                         function(x) unlist(lapply(x@Polygons, function(y)
                                                   .arcp(y@coords)))))
    lhol <- unlist(lapply(polygons(spdf)@polygons,
                          function(x) unlist(lapply(x@Polygons, function(y)
                                                    y@hole))))
    sum(lar[!lhol])-sum(lar[lhol])
}
