"lonlat2xy" <-
function (lnlt, miles = FALSE) {
  # converts lon/lat coordinates to x/y coordinates, using a very simple projection
  # the projection calculates, for all points, the great circle distance in the x direction to the mean longitude and in the y direction to the mean latitude, and uses these distances as the x-y coordinate of the location
  # require('fields')  # no longer needed as I added rdist.earth to spectralGP
  if(!(is.matrix(lnlt) || is.data.frame(lnlt)) || nrow(lnlt)<1 || ncol(lnlt)!=2){
    stop("lnlt must be a matrix or dataframe with two columns, where the first column indicates the x-coordinate and the second the y-coordinate")
  }
  if(min(lnlt[,1])<(-360) || max(lnlt[,1])>360 || min(lnlt[,2])<(-90) || max(lnlt[,1])>90){
    stop("lnlt values must lie in (-360,360)X(-90,90)")
  }
  n <- nrow(lnlt)
  ctr <- c(mean(lnlt[, 1]), mean(lnlt[, 2]))
  ctx <- cbind(lnlt[, 1], rep(ctr[2], n))
  cty <- cbind(rep(ctr[1], n), lnlt[, 2])
  dstx <- c(rdist.earth(matrix(ctr, 1, 2), ctx, miles = FALSE))
  dsty <- c(rdist.earth(matrix(ctr, 1, 2), cty, miles = FALSE))
  cbind(dstx * (2 * (lnlt[, 1] > ctr[1]) - 1), dsty * (2 * (lnlt[, 2] > ctr[2]) - 1))
}
