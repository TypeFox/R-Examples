# This is a hidden function of the l2boost package.

# plots.lines is used by \code{\link{plot.l2boost}} to the path lines (each j, against each r-step)
# 
# @param xval vector of x-values corresponding to the path y-values (default: NULL index of path)
# @param ind Coordinate of the path (for coloring individual paths)
# @param path Plot the path values along the y-axis
# @param l.crit change the color at the value of m=l.crit
# @param active active set coloring (default: TRUE)
# @param col vector of color values lnegth >= 1 (default: NULL use built in scheme)
# 

plot.lines <- function (xval = NULL, ind, path, l.crit, active = TRUE, col=NULL) {
  M <- length(path)
  
  actCount = 1
  lty=2
  if (is.null(xval)) xval <- 1:M
  for (k in ind) {
    path.k <- sapply(1:M, function(m){path[[m]][k]})
    if (active) {    
      when.k <- which(l.crit == k)
      if(is.null(col)){
        active.k <- sapply(1:M, function(m) {
          active.M <- 1
          if ((length(when.k) == 0) | (length(when.k) > 0 && min(when.k) > m)) {
            active.M <- 0
          }
          active.M})
        #     points(xval, path.k, col = col.k, pch = 16, cex = 0.5)
        path.mod.k <- path.k
        path.mod.k[active.k == 0] <- NA
        lines(xval, path.mod.k, lty = 2, col = "red") 
        path.mod.k <- path.k
        path.mod.k[active.k != 0] <- NA
        lines(xval, path.mod.k, lty = 2, col = "blue") 
      }else{
        if(length(col) > 1){
          if(actCount > length(col) ) actCount = 1
          clr = col[actCount]
          actCount = actCount+1
          lty=1
        }else{
          clr= col
          lty=2
        }
        lines(xval, path.k, lty = lty, lwd = 3,col = clr)
      }
    }else {
      lines(xval, path.k, lty = 1, col = "gray")
    }
  }
}
