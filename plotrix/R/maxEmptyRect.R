# function to find the largest rectangle not containing any of the points
# specified by x and y
# adapted by Hans Borchers from 
# A. Naamad, D. T. Lee, and W.-L. Hsu (1984). On the Maximum Empty
# Rectangle Problem. Discrete Applied Mathematics, Vol. 8, pp. 267--277.

maxEmptyRect <- function(ax, ay, x, y) {
 n <- length(x)
 d <- sort(c(ax, x))
 D <- diff(d)
 m <- which.max(D)
 # check vertical slices
 mgap <- D[m]
 maxr <- mgap * (ay[2] - ay[1])
 maxR <- c(d[m], ay[1], d[m+1], ay[2])
 o <- order(y)
 X <- x[o]; Y <- y[o]
 for (i in 1:n) {
  tl <- ax[1]; tr <- ax[2]
  if (i < n) {
   for (j in (i+1):n) {
    if (X[j] > tl && X[j] < tr) {
     # check horizontal slices (j == i+1)
     # and (all) rectangles above (X[i], Y[i])
     area <- (tr-tl)*(Y[j]-Y[i])
     if (area > maxr) {
      maxr <- area
      maxR <- c(tl, Y[i], tr, Y[j])
     }
     if (X[j] > X[i]) tr <- X[j]
     else tl <- X[j]
    }
   }
  }
  # check open rectangles above (X[i], Y[i])
  area <- (tr-tl)*(ay[2]-Y[i])
  if (area > maxr) {
   maxr <- area
   maxR <- c(tl, Y[i], tr, ay[2])
  }
 }
 for (i in 1:n) {
  # check open rectangles above (X[i], Y[i])
  ri <- min(ax[2], X[Y > Y[i] & X > X[i]])
  li <- max(ax[1], X[Y > Y[i] & X < X[i]])
  area <- (ri-li)*(ay[2]-Y[i])
  if (area > maxr) {
   maxr <- area
   maxR <- c(li, Y[i], ri, ay[2])
  }
  # check open rectangles below (X[i], Y[i])
  ri <- min(ax[2], X[Y < Y[i] & X > X[i]])
  li <- max(ax[1], X[Y < Y[i] & X < X[i]])
  area <- (ri-li)*(Y[i]-ay[1])
  if (area > maxr) {
   maxr <- area
   maxR <- c(li, ay[1], ri, Y[i])
  }
 }
 xcenter<-sum(maxR[c(1,3)])/2
 ycenter<-sum(maxR[c(2,4)])/2
 return(list(area = maxr, rect = maxR, x = xcenter, y = ycenter))
}
