library(rgdal)
library(spdep)
library(rgeos)


################################################## 

# GEOGRAPHIC FUNCTIONS

################################################## 


#' A function to compute the spatial units' areas
#'
#' @usage area(spatobj = NULL, folder = NULL, shape = NULL)
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a Spatial object (SpatialPolygonsDataFrame)
#' @return A vector with spatial units' areas. 
#' @description The function is based on \pkg{rgdal} package and 
#' can be used by providing a shape file or a R spatial object 
#' (SpatialPolygonsDataFrame).
#' @examples  area(GreHSize) 
#' 
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' area(folder = foldername, shape = shapename)
#' @seealso  Other spatial functions used for segregation indices 
#' computation: \code{\link{contig}}, \code{\link{perimeter}}, 
#' \code{\link{distance}}, \code{\link{distcenter}}, 
#' \code{\link{boundaries}}, \code{\link{xgeo}}
#' @export

area <- function(spatobj = NULL, folder = NULL, shape = NULL) {
    if (is.null(spatobj)) 
        spatobj <- rgdal::readOGR(dsn = folder, layer = shape)
    provi <- slot(spatobj, "polygons")
    area <- sapply(provi, slot, "area")
    return(area)
}

#' A function to compute the contiguity matrix
#'
#' @usage contig(spatobj = NULL, folder = NULL, shape = NULL)
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a Spatial object (SpatialPolygonsDataFrame)
#' @return A first order contiguity (adjacency) matrix, where each 
#' element [\emph{i,j}] equals 1 if \emph{i}-th  and \emph{j}-th  
#' spatial units are adjacent, 0 otherwise (queen criteria)
#' @description The function is based on \pkg{rgdal} and 
#' \pkg{spdep} packages and it can be used by providing a shape 
#' file or a R spatial object (SpatialPolygonsDataFrame).
#' @examples  contig(GreHSize) 
#' 
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' contig(folder = foldername, shape = shapename)
#' @seealso  Other spatial functions used for segregation indices 
#' computation: \code{\link{area}}, \code{\link{perimeter}}, 
#' \code{\link{distance}}, \code{\link{distcenter}}, 
#' \code{\link{boundaries}}, \code{\link{xgeo}}
#' @export

contig <- function(spatobj = NULL, folder = NULL, shape = NULL) {
    if (is.null(spatobj)) 
        spatobj <- rgdal::readOGR(dsn = folder, layer = shape)
    data_ngb <- spdep::poly2nb(spatobj)
    contig <- spdep::nb2mat(data_ngb, style = "B")
    return(contig)
}


#' A function to compute the spatial units' perimeters
#'
#' @usage perimeter(spatobj = NULL, folder = NULL, shape = NULL)
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a Spatial object (SpatialPolygonsDataFrame)
#' @return A vector with spatial units' perimeters. 
#' @description The function is based on \pkg{rgdal} and \pkg{rgeos} 
#' packages and it can be used by providing a shape 
#' file or a R spatial object (SpatialPolygonsDataFrame).
#' @examples  perimeter(GreHSize) 
#' 
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' perimeter(folder = foldername, shape = shapename)
#' @seealso  Other spatial functions used for segregation indices 
#' computation:  \code{\link{area}}, \code{\link{contig}}, 
#' \code{\link{distance}}, \code{\link{distcenter}}, 
#' \code{\link{boundaries}}, \code{\link{xgeo}}
#' @export

perimeter <- function(spatobj = NULL, folder = NULL, shape = NULL) {
    if (is.null(spatobj)) 
        spatobj <- rgdal::readOGR(dsn = folder, layer = shape)
    perim <- vector(length = nrow(spatobj))
    for (i in 1:nrow(spatobj)) 
      perim[i] <- rgeos::gLength(spatobj[i, ])
    return(perim)
}


#' A function to compute the distance matrix between centroids 
#' of spatial units
#'
#' @usage distance(spatobj = NULL, folder = NULL, shape = NULL)
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a Spatial object (SpatialPolygonsDataFrame)
#' @return A matrix with the distance between spatial units centroids
#' @description The function is based on \pkg{rgdal} and \pkg{rgeos} 
#' packages and it can be used by providing a  shape file or a R 
#' spatial object (SpatialPolygonsDataFrame).
#' @examples  distance(GreHSize) 
#' 
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' distance(folder = foldername, shape = shapename)
#' @seealso  Other spatial functions used for segregation indices 
#' computation: \code{\link{area}}, \code{\link{contig}}, 
#' \code{\link{perimeter}}, \code{\link{distcenter}}, 
#' \code{\link{boundaries}}, \code{\link{xgeo}}
#' @export

distance <- function(spatobj = NULL, folder = NULL, shape = NULL) {
    if (is.null(spatobj)) 
        spatobj <- rgdal::readOGR(dsn = folder, layer = shape)
    dist <- matrix(0, nrow = nrow(spatobj), ncol = nrow(spatobj))
    centroids <- vector("list", nrow(spatobj))
    for (i in 1:nrow(spatobj)) centroids[[i]] <- rgeos::gCentroid(spatobj[i, ])
    for (i in 1:(nrow(spatobj) - 1)) 
      for (j in (i + 1):nrow(spatobj)) 
        dist[i, j] <- rgeos::gDistance(centroids[[i]], 
        centroids[[j]])
    dist <- dist + t(dist)
    return(dist)
}

#' A function to compute the distance from spatial units centroids 
#' to the center
#'
#' @usage distcenter(spatobj = NULL, folder = NULL, shape = NULL, center = 1)
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a Spatial object (SpatialPolygonsDataFrame)
#' @param center - the row number of the center
#' @return A vector with the distance to the center's centroid
#' @description The function is based on \pkg{rgdal} and \pkg{rgeos} 
#' packages and it can be used by providing a shape file or a R 
#' spatial object (SpatialPolygonsDataFrame).
#' @examples  distcenter(GreHSize, center = 19) 
#' 
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' distcenter(folder = foldername, shape = shapename, center = 19)
#' @seealso  Other spatial functions used for segregation indices 
#' computation: \code{\link{area}}, \code{\link{contig}}, 
#' \code{\link{perimeter}}, \code{\link{distance}}, 
#' \code{\link{boundaries}}, \code{\link{xgeo}}
#' @export


distcenter <- function(spatobj = NULL, folder = NULL, 
                       shape = NULL, center = 1) {
    if (is.null(spatobj)) 
        spatobj <- rgdal::readOGR(dsn = folder, layer = shape)
    distcenter <- vector(length = nrow(spatobj))
    centroids <- vector("list", nrow(spatobj))
    for (i in 1:nrow(spatobj)) 
      centroids[[i]] <- rgeos::gCentroid(spatobj[i, ])
    for (i in 1:nrow(spatobj)) 
      distcenter[i] <- rgeos::gDistance(centroids[[i]], centroids[[center]])
    return(distcenter)
}

#' A function to compute the matrix of common boundaries
#'
#' @usage boundaries(spatobj = NULL, folder = NULL, shape = NULL)
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a Spatial object (SpatialPolygonsDataFrame)
#' @return A matrix with common boundaries between spatial units
#' @description The function is based on \pkg{rgdal} and \pkg{rgeos} 
#' packages and it can be used by providing a shape file 
#' or a R spatial object (SpatialPolygonsDataFrame).
#' @examples  boundaries(AnnHAge) 
#' 
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'AnnHAge'
#' boundaries(folder = foldername, shape = shapename)
#' @seealso  Other spatial functions used for segregation indices 
#' computation:  \code{\link{area}}, \code{\link{contig}}, 
#' \code{\link{perimeter}}, \code{\link{distance}}, 
#' \code{\link{distcenter}}, \code{\link{xgeo}}
#' @export


boundaries <- function(spatobj = NULL, folder = NULL, shape = NULL) {
    if (is.null(spatobj)) 
        spatobj <- rgdal::readOGR(dsn = folder, layer = shape)
    boundaries <- matrix(0, nrow = nrow(spatobj), ncol = nrow(spatobj))
    for (i in 1:(nrow(spatobj) - 1)) for (j in (i + 1):nrow(spatobj)) {
        provi <- rgeos::gIntersection(spatobj[i, ], spatobj[j, ])
        if (!is.null(provi)) 
            boundaries[i, j] <- rgeos::gLength(provi)
    }
    boundaries <- boundaries + t(boundaries)
    return(boundaries)
}

#' A function to compute all spatial info for segregation indexes
#'
#' @usage xgeo(spatobj = NULL, folder = NULL, shape = NULL, center = 1)
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a Spatial object (SpatialPolygonsDataFrame)
#' @param center - the row number of the center
#' @return A list that contains all geographic information 
#' needed to calculate segregation indexes.
#' @description The function is based on \pkg{rgdal}, \pkg{rgeos} 
#' and \pkg{spdep} packages and it can be used by providing a 
#' shape file or a R spatial object (SpatialPolygonsDataFrame).
#' @examples  xgeo(AnnHAge, center = 1) 
#' 
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'AnnHAge'
#' xgeo(folder = foldername, shape = shapename, center = 19)
#' @seealso  Other spatial functions used for segregation indices 
#' computation:  \code{\link{area}}, \code{\link{contig}}, 
#' \code{\link{perimeter}}, \code{\link{distance}}, 
#' \code{\link{distcenter}}, \code{\link{boundaries}}
#' @export

xgeo <- function(spatobj = NULL, folder = NULL, shape = NULL, center = 1) {
    if (is.null(spatobj)) 
        spatobj <- rgdal::readOGR(dsn = folder, layer = shape)
    data_ngb <- spdep::poly2nb(spatobj)
    contig <- spdep::nb2mat(data_ngb, style = "B")
    provi <- slot(spatobj, "polygons")
    area <- sapply(provi, slot, "area")
    perim <- vector(length = nrow(spatobj))
    for (i in 1:nrow(spatobj)) 
      perim[i] <- rgeos::gLength(spatobj[i, ])
    dist <- matrix(0, nrow = nrow(spatobj), ncol = nrow(spatobj))
    centroids <- vector("list", nrow(spatobj))
    for (i in 1:nrow(spatobj)) 
      centroids[[i]] <- rgeos::gCentroid(spatobj[i, ])
    for (i in 1:(nrow(spatobj) - 1)) 
      for (j in (i + 1):nrow(spatobj)) 
        dist[i, j] <- rgeos::gDistance(centroids[[i]], 
        centroids[[j]])
    dist <- dist + t(dist)
    distcenter <- vector(length = nrow(spatobj))
    for (i in 1:nrow(spatobj)) 
      distcenter[i] <- rgeos::gDistance(centroids[[i]], centroids[[center]])
    boundaries <- matrix(0, nrow = nrow(spatobj), ncol = nrow(spatobj))
    for (i in 1:(nrow(spatobj) - 1)) 
      for (j in (i + 1):nrow(spatobj)) {
        provi <- rgeos::gIntersection(spatobj[i, ], spatobj[j, ])
        if (!is.null(provi)) 
            boundaries[i, j] <- rgeos::gLength(provi)
    }
    boundaries <- boundaries + t(boundaries)
    resultlist <- list(perimeter = perim, area = area, contiguity = contig, 
          boundaries = boundaries, distance = dist, distcenter = distcenter)
    return(resultlist)
} 
