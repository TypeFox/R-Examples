"rdist.earth" <-
function (loc1, loc2, miles = TRUE, R = NULL) {
  # this is the fields function rdist.earth()
  # calculates great circle distances between all pairs of points
    if (is.null(R)) {
        if (miles) 
            R <- 3963.34
        else R <- 6378.388
    }
    if (missing(loc2)) 
        loc2 <- loc1
    coslat1 <- cos((loc1[, 2] * pi)/180)
    sinlat1 <- sin((loc1[, 2] * pi)/180)
    coslon1 <- cos((loc1[, 1] * pi)/180)
    sinlon1 <- sin((loc1[, 1] * pi)/180)
    coslat2 <- cos((loc2[, 2] * pi)/180)
    sinlat2 <- sin((loc2[, 2] * pi)/180)
    coslon2 <- cos((loc2[, 1] * pi)/180)
    sinlon2 <- sin((loc2[, 1] * pi)/180)
    pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*% 
        t(cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2))
    R * acos(ifelse(pp > 1, 1, pp))
}
