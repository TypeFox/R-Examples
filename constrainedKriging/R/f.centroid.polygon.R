f.centroid.polygon <- function(t.polygon)

### purpose: calculate the  centroid of polygons
###          function gives a list back each element is a
###          1 X 2 matrix with the coords of the centroid
###
### arguments:
###            t.polygon = polygon list Class: "SpatialPolygons" "Spatial"   
###
### author: ch.hofer       
### date: 20.2.2006
{
# as.matrix.polygon and centroid.polygon are internally function of the maps packag
as.matrix.polygon <- function(x, ...) {
  p = x
  if(is.null(p)) return(p)
  if(is.list(p) && !is.data.frame(p)) p <- cbind(p$x, p$y)
  p
}

centroid.polygon <- function(p) {
  if(is.null(p)) return(c(NA, NA))
  p <- as.matrix.polygon(p)
  n <- nrow(p)
  x1 <- p[, 1]
  i2 <- c(n, 1:(n - 1))
  x2 <- p[i2, 1]
  y1 <- p[, 2]
  y2 <- p[i2, 2]
  a <- x1*y2 - x2*y1
  s <- sum(a)*3
  if(s == 0) c(mean(x1), mean(y1))
  else c(sum((x1 + x2)*a)/s, sum((y1 + y2)*a)/s)
}

    
t.centroids <-lapply(t.polygon@polygons,
                          function(x){
                            return(
                                   centroid.polygon(x@Polygons[[1]]@coords)
                                   ) ## end return
                          } ## end function
                     ) ## end lapply

### change list into a matrix dim n x 2 where n is the number of polygons
t.x <- rep(c(T,F), length(t.centroids))
t.centroids <- cbind(x.centroid  = unlist(t.centroids)[t.x == T],
                     y.centroid = unlist(t.centroids)[t.x == F])
          
return(t.centroids)
}
