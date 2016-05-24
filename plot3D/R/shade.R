## =============================================================================
## =============================================================================
## Shaded colors for 3-D images
## =============================================================================
## =============================================================================

## =============================================================================
## Calculate surface normals from (x, y, z) matrices
## =============================================================================

normal.matrix <- function(x, y, z, Extend = TRUE, na.rm = FALSE) {  # the x- y- and z matrices
 
 # x, y and z: matrices of same dimension
  if (Extend & !na.rm) {
    xx <- extend(x) 
    yy <- extend(y) 
    zz <- extend(z) 
  } else if (Extend & na.rm) {
    xx <- extend.na(x) 
    yy <- extend.na(y) 
    zz <- extend.na(z) 
  } 
  else {
    xx <- x 
    yy <- y 
    zz <- z 
  }

 # the facets:
  ii <- 1 : (nrow(xx)-1)
  jj <- 1 : (ncol(xx)-1)
    
 # normals in x, y, and z direction are matrices
  N.x <- N.y <- N.z <- matrix(nrow = nrow(xx)-1, ncol = ncol(xx)-1, data = NA)
    
 # choose points on each vertex
  for (i in ii) {
    ip1 <- cbind(i,     jj+1)
    ip2 <- cbind(i+1  , jj)
    ip3 <- cbind(i+1  , jj+1)
    ip4 <- cbind(i    , jj)

    p1  <- cbind(xx[ip1] , yy[ip1], zz[ip1])
    p2  <- cbind(xx[ip2] , yy[ip2], zz[ip2])
    p3  <- cbind(xx[ip3] , yy[ip3], zz[ip3])
    p4  <- cbind(xx[ip4] , yy[ip4], zz[ip4])

   # two vectors to represent these points
    V1 <- p2 - p1
    V2 <- p3 - p4

   # the (unnormalised) normals
    N.x [i, ] <-  V1[, 2]*V2[, 3] - V1[, 3]*V2[, 2]
    N.y [i, ] <- -V1[, 1]*V2[, 3] + V1[, 3]*V2[, 1]
    N.z [i, ] <-  V1[, 1]*V2[, 2] - V1[, 2]*V2[, 1]
  }
  # normalise
  Norm <- sqrt(N.x^2 + N.y^2 + N.z^2) # normalisation factor
  Norm[Norm == 0] <-1
  list (u = N.x/Norm, v = N.y/Norm, w = N.z/Norm)  
}

## =============================================================================

normal.points <- function(p1, p2, p3, p4) {  #x, y, z of 4 pts
   # two vectors to represent these points
  V1 <- p2 - p1
  V2 <- p3 - p4

   # the (unnormalised) normals
  N.x  <-  V1[2, ]*V2[3, ] - V1[3, ]*V2[2, ]
  N.y  <- -V1[1, ]*V2[3, ] + V1[3, ]*V2[1, ]
  N.z  <-  V1[1, ]*V2[2, ] - V1[2, ]*V2[1, ]

  # normalise
  Norm <- sqrt(N.x^2 + N.y^2 + N.z^2) # normalisation factor
  Norm[Norm == 0] <-1
  list (u = N.x/Norm, v = N.y/Norm, w = N.z/Norm)  
}

## =============================================================================

normal.points.tri <- function(p1, p2, p3) {  #x, y, z of 3 pts
   # two vectors to represent these points
  V1 <- p2 - p1
  V2 <- p3 - p1

   # the (unnormalised) normals
  N.x  <-  V1[2, ]*V2[3, ] - V1[3, ]*V2[2, ]
  N.y  <- -V1[1, ]*V2[3, ] + V1[3, ]*V2[1, ]
  N.z  <-  V1[1, ]*V2[2, ] - V1[2, ]*V2[1, ]

  # normalise
  Norm <- sqrt(N.x^2 + N.y^2 + N.z^2) # normalisation factor
  Norm[Norm == 0] <-1
  list (u = N.x/Norm, v = N.y/Norm, w = N.z/Norm)  
}

## =============================================================================
## Setup light based on light angles
## =============================================================================
# A translation to R from the C-code in plot3d.c 

setuplight <- function(phil, thetal) {

 # rotation in x-direction
  cosp <- cos(-phil/180*pi)
  sinp <- sin(-phil/180*pi)
  rotX <- matrix(nrow = 4, data = c(1, 0,   0,      0, 
                                    0, cosp, -sinp, 0, 
                                    0, sinp,  cosp, 0,
                                    0, 0,   0,      1))
  VT <- diag(nrow = 4)  %*% rotX

 # rotation in z-direction 
  cost <- cos(thetal/180*pi)
  sint <- sin(thetal/180*pi)
  rotZ <- matrix(nrow = 4, data = c(cost, -sint, 0, 0,
                                    sint,  cost, 0, 0,
                                    0,   0,      1, 0,
                                    0,   0,      0, 1))
  VT <- VT %*% rotZ

  light <-  c(0, -1, 0, 1) %*% VT 
  return(light)

}

## =============================================================================
## Create 3-D facet colors with shading or lighting. 
## =============================================================================

facetcols.tri <- function(tri, col, shade){

  Nr <- nrow(tri)/3
  A <- array(dim = c(3, 3, Nr), data = t(tri))
  
  if (length(col) != Nr)
    col <- rep(col, length.out = Nr)
  A[1,,] <-  A[1,,] *shade$xs
  A[2,,] <-  A[2,,] *shade$ys
  A[3,,] <-  A[3,,] *shade$zs
  light   <- setuplight(shade$lphi, shade$ltheta) [1:3]
  Normals <- normal.points.tri(A[,1,], A[,2,], A[,3,])
    
  return(facetcols.shadelight(light, Normals, col, shade))
}

## =============================================================================
## x-y-z is a matrix
## =============================================================================

facetcols <- function(x, y, z, col, shade, Extend = TRUE){

 # +90 for "rotation to horizontal"
  if (Extend & length(col) != length(x))
    col <- rep(col, length.out = length(x))

  else if (!Extend & length(col) != prod(dim(x)-1))
    col <- rep(col, length.out = prod(dim(x)-1))

  light   <- setuplight(shade$lphi, shade$ltheta) [1:3]
  Normals <- normal.matrix(x * shade$xs, y * shade$ys, z * shade$zs, Extend)
    
  return(facetcols.shadelight(light, Normals, col, shade))
}

## =============================================================================
## facet colors with transparancy
## =============================================================================

facetcols.shadelight <- function(light, Normals, col, shade){
  
  # we keep "transparent" colors
  ii <- which (col == "transparent")
   
  if (shade$type == "shade")
   Col <- facetcols.shade(light, Normals, col, shade$shade)

  else if (shade$type == "light")
   Col <- facetcols.light(light, Normals, col, shade)
  
  if (! is.null(shade$alpha))
    Col <- setalpha(Col, shade$alpha)

  if (length(ii) > 0)
    Col[ii] <- "transparent"
      
  return(Col)  
}

## =============================================================================
## facet colors with simplified phong lighting. 
## =============================================================================

facetcols.light <- function(light, Normals, col, shade) {

 # defaults
  p <- list(ambient = 0.3, diffuse = 0.6, specular = 1.,
            exponent = 20, sr = 0, alpha = 1)
  nmsC <- names(p)
  p[(namc <- names(shade))] <- shade

 # this is different from shaded colors - 
#  Sum <- Normals$u*light[1] + Normals$v*light[2] + Normals$w*light[3]
 # use same as shaded.colors
  Sum <- 0.5*(Normals$u*light[1] + Normals$v*light[2] + Normals$w*light[3] +1)

  Is <- as.vector(p$specular * abs(Sum) ^ p$exponent)
  Id <- as.vector(p$diffuse  * pmax(Sum, 0))

  rgbcol  <- t(col2rgb(col) / 255)
  Lrgbcol <- pmin((p$ambient + Id + p$sr * Is) * rgbcol + (1 - p$sr) * Is, 1)

  Lrgbcol[is.na(Lrgbcol)] <- 0
  if (is.null(p$alpha))
    p$alpha <- 1  # necessary for R < 3.0
  rgb(Lrgbcol[,1], Lrgbcol[,2], Lrgbcol[,3], p$alpha)

}

## =============================================================================
## facet colors with shading
## =============================================================================

facetcols.shade <- function(light, Normals, col, shade){
  
  if (is.na(shade))
    return(col)
  shade <- abs(shade) 
  if (shade < 0 | shade > 1) 
    stop("'shade' should be a value inbetween 0 and 1")
 
  Sum <- 0.5*(Normals$u*light[1] + Normals$v*light[2] + Normals$w*light[3] +1)
  Shade <- Sum^shade
  Shade[is.na(Shade)] <- 0
  
  RGB <- t(col2rgb(col)) * as.vector(Shade) / 255
  col[] <- rgb(RGB)      # alpha = 1

  return(col)
}


facetcolsImage <- function (x, y, z, xlim, ylim, zlim, shade, lighting, alpha, 
  ltheta, lphi, Col, NAcol) {

  if (is.null(xlim))
    xlim <- range(x, na.rm = TRUE)
  if (is.null(ylim))
    ylim <- range(y, na.rm = TRUE)
  if (is.null(zlim))
    zlim <- range(z, na.rm = TRUE)
  xs <- 0.5 *abs(diff(xlim))
  ys <- 0.5 *abs(diff(ylim))
  zs <- 0.5 *abs(diff(zlim))
  xs <- ifelse (xs == 0, 1, 1 / xs)
  ys <- ifelse (ys == 0, 1, 1 / ys)
  zs <- ifelse (zs == 0, 1, 1 / zs)

  if (! is.matrix(x)) {
    xy <- mesh(x,y)
    x <- xy$x
    y <- xy$y  
  }
  
  light   <- setuplight(lphi, ltheta) [1:3]
  Normals <- normal.matrix(x * xs, y * ys, z * zs, Extend = TRUE, 
    na.rm = any(is.na(z)))
  ina <- which (is.na(Normals$u))
    
  List <- list (shade = shade, alpha = alpha, type = "none")

  if (! is.null(lighting)) {         
    if (is.character(lighting))
      List$type <- "light"
    else if (is.logical(lighting)) {
      if (lighting)
        List$type <- "light"
    } else if (is.list(lighting)) {
      if (!is.null(lighting$type)) 
        List$type <- lighting$type
      else   
        List$type <- "light"
      lighting$type <- NULL
      List <- c(List, lighting)
    }
  }
  if (! is.null(shade))
    if (! is.na(shade) & List$type == "none") # lighting overrules shade
      List$type <- "shade"     
  if (is.null(List$shade)) 
    List$shade <- NA
  
  # we keep "transparent" colors
  col <- Col
  ii <- which (col == "transparent")
   
  Col[] <- facetcols.shadelight(light, Normals, col, List)

  if (! is.null(alpha))
    Col <- setalpha(Col, alpha)

  if (length(ina) > 0)
    Col[ina] <- NAcol
    
  if (length(ii) > 0 )
    Col[ii] <- "transparent"

  return(Col)

}
