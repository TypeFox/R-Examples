
"distSpeed" <- function(pt1, pt2, method=c("Meeus", "VincentyEllipsoid"))
{
    ## Value: A 3-column matrix with distance, time elapsed and speed
    ## between two points or set of points.
    ## --------------------------------------------------------------------
    ## Arguments: pt1 and pt2=matrices for each point, with three columns;
    ## the first for a POSIXct object with time for each point, the second
    ## for longitude, and the third for latitude.  method=character; which
    ## of the distance algorithms from geosphere package to use (only
    ## default parameters used).
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (nrow(pt1) != nrow(pt2)) {
        stop("pt1 and pt2 must have the same number of rows")
    } else if (ncol(pt1) != 3 || ncol(pt2) != 3) {
        stop("pt1 and pt2 must both have 3 columns")
    } else if (nrow(pt1) < 1 || ncol(pt2) < 1) {
        stop("pt1 and pt2 must each have at least 1 row")
    }
    method <- match.arg(method)
    switch(method,
           Meeus = {
               distance <- geosphere::distMeeus(pt1[, 2:3], pt2[, 2:3])
           },
           VincentyEllipsoid = {
               distance <- geosphere::distVincentyEllipsoid(pt1[, 2:3],
                                                            pt2[, 2:3])
           })
    pt1[, 1] <- as.numeric(pt1[, 1])
    pt2[, 1] <- as.numeric(pt2[, 1])
    ## Distance (in Km)
    distance <- distance / 1000
    ## Calculate time difference (in seconds) between locations.
    timdiff <- abs(pt1[, 1] - pt2[, 1])
    ## Speed in m/s.
    speed <- ifelse(timdiff == 0, 0, (distance * 1000) / timdiff)
    cbind(distance, time.elapsed=timdiff, speed)
}
