
################################################## 

# SEGREGATION INDEXES

################################################## 

# EVENESS INDEXES

################### 


#' A function to compute Duncan & Duncan segregation index
#'
#' @usage Duncan (x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A vector with Duncan segregation index 
#' @references Duncan O. D. and Duncan B. (1955) \emph{A 
#' Methodological Analysis of Segregation Indexes}. 
#' American Sociological Review 41, pp. 210-217
#' @description Duncan's segregation index is one of the most used 
#' seggregation index and varies from 0 (perfect homogenous spatial 
#' distribution of a group) to 1 (maximum spatial segregation)
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' Duncan(x) 
#' @seealso  Other indices based on Duncan's index: 
#' \code{\link{Wong}}, \code{\link{Morill}}
#' @seealso Other evenness intragroup  indices: 
#' \code{\link{Gini}}, \code{\link{Gorard}}
#' @seealso Intergroup dissimilarity index: 
#' \code{\link{DI}}
#' @export


Duncan <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    tx <- rowSums(x)
    px <- x/tx
    for (i in 1:ncol(x)) {
        px[, i] <- tx * abs(px[, i] - pTotal[i])
        result[i] <- sum(px[, i])/(2 * Total * pTotal[i] * (1 - pTotal[i]))
    }
    return(result)
}


#' A function to compute Morill's segregation index 
#'
#' @usage Morill(x, c = NULL, spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param c - a binary contiguity (adjacency) symetric matrix where each 
#' element \emph{Cij} equals 1 if \emph{i}-th and \emph{j}-th spatial units 
#' are adjacent, 0 otherwise. 
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A vector with Morill segregation index 
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' contiguity <- contig(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' Morill(x, c = contiguity) 
#' 
#' Morill(x, spatobj = GreHSize)
#' 
#' Morill(x, folder = foldername, shape = shapename) 
#' 
#' @references Morill B. (1991) \emph{On the measure of geographic 
#' segregation}. Geography research forum, 11, pp. 25-36.
#' @description Morill's segregation index develops \code{\link{Duncan}}'s 
#' index by taking into account the interactions between spatial units 
#' (contiguity). The function can be used in two ways: by providing a 
#' contiguity matrix or a geographic source (spatial object or shape file) 
#' wich will be used to compute the contiguity matrix within the function
#' @seealso Other evenness intragroup  indices: 
#' \code{\link{Duncan}}, \code{\link{Wong}}, \code{\link{Gini}}, 
#' \code{\link{Gorard}} 
#' @seealso Intergroup dissimilarity index: \code{\link{DI}}
#' @export


Morill <- function(x, c = NULL, spatobj = NULL, folder = NULL, 
                   shape = NULL) {
    x <- as.matrix(x)
    if (is.null(c)) 
        c <- contig(spatobj = spatobj, folder = folder, shape = shape)
    IS <- vector(length = ncol(x))
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    tx <- rowSums(x)
    px <- x/tx
    for (i in 1:ncol(x)) {
        px[, i] <- tx * abs(px[, i] - pTotal[i])
        IS[i] <- sum(px[, i])/(2 * Total * pTotal[i] * (1 - pTotal[i]))
    }
    result <- vector(length = ncol(x))
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    t <- rowSums(x)
    p <- x/t
    for (k in 1:ncol(x)) {
        for (i in 1:nrow(x)) pij[k, , i] <- p[i, k]
        for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - p[i, k])
        matprovi <- c %*% pij[k, , ]
        matprovi <- matprovi * diag(nrow(x))
        result[k] <- IS[k] - sum(matprovi)/sum(c)
    }
    return(result)
}


#' A function to compute Wong's segregation index 
#'
#' @usage Wong(x, b = NULL, p = NULL, a = NULL,  spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param a - a vector with areas of spatial units
#' @param b - a common boundaries matrix where each element \emph{Bij} 
#' equals the shared boundary of \emph{i}-th and \emph{j}-th spatial 
#' units. 
#' @param p - a vector with the perimeters of spatial units
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A vector with Wong segregation index 
#' @examples x <- slot(AnnHAge, 'data')[ ,3:5]
#' bound <- boundaries(AnnHAge)
#' per <- perimeter(AnnHAge)
#' ar <- area(AnnHAge)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'AnnHAge'
#' 
#' Wong(x, b = bound, p = per, a = ar) 
#' 
#' Wong(x, spatobj = AnnHAge)
#' 
#' Wong(x, folder = foldername, shape = shapename) 
#' 
#' @references Wong D. W. S. (1998) \emph{Measuring multiethnic spatial 
#' segregation}. Urban Geography, 19 (1), pp. 77-87.
#' @description Wong's segregation index develops \code{\link{Duncan}}'s 
#' index by taking into account the interaction between spatial units 
#' (common boundaries and perimeters). The function can be used 
#' in two ways: by providing the boundaries matrix and perimeter 
#' vector or a geographic source (spatial object or shape file) wich 
#' will be used to compute the geographic information within the function
#' @seealso Other evenness intragroup  indices: \code{\link{Duncan}}, 
#' \code{\link{Morill}}, \code{\link{Gini}}, \code{\link{Gorard}}
#' @seealso Intergroup dissimilarity index: \code{\link{DI}}
#' @export


Wong <- function(x, b = NULL, p = NULL, a = NULL, spatobj = NULL, 
                 folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(b)) 
        b <- boundaries(spatobj = spatobj, folder = folder, shape = shape)
    if (is.null(p)) 
        p <- perimeter(spatobj = spatobj, folder = folder, shape = shape)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    IS <- vector(length = ncol(x))
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    tx <- rowSums(x)
    px <- x/tx
    for (i in 1:ncol(x)) {
        px[, i] <- tx * abs(px[, i] - pTotal[i])
        IS[i] <- sum(px[, i])/(2 * Total * pTotal[i] * (1 - pTotal[i]))
    }
    result <- vector(length = ncol(x))
    PerAij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    t <- rowSums(x)
    pp <- x/t
    maxPerA <- max(p/a)
    w <- b/p
    for (k in 1:ncol(x)) {
        for (i in 1:nrow(x)) pij[k, , i] <- pp[i, k]
        for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - pp[i, k])
        for (i in 1:nrow(x)) PerAij[k, , i] <- p[i]/a[i]
        for (i in 1:nrow(x)) PerAij[k, i, ] <- PerAij[k, i, ] + p[i]/a[i]
        matprovi <- (w/sum(w)) %*% (pij[k, , ] * PerAij[k, , ])
        matprovi <- matprovi * diag(nrow(x))
        result[k] <- IS[k] - 0.25 * sum(matprovi)/maxPerA
    }
    return(result)
}


#' A function to compute Gorard's segregation index 
#'
#' @usage Gorard(x)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A vector with Gorard segregation index 
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' Gorard(x)
#' @references Gorard S. (2000) \emph{Education and Social Justice}. 
#' Cardiff, University of Wales Press
#' @description Gorard's index is an alternative to segregation indices 
#' based on \code{\link{Duncan}}'s index. The index varies between 0 
#' (minimum segregation) and 1 (maximum segregation).
#' @seealso Other evenness intragroup  indices: \code{\link{Duncan}}, 
#' \code{\link{Morill}}, \code{\link{Wong}}, \code{\link{Gini}}
#' @seealso Intergroup dissimilarity index: \code{\link{DI}}
#' @export


Gorard <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    tx <- rowSums(x)
    varTotal <- colSums(x)
    Total <- sum(x)
    for (k in 1:ncol(x)) 
      result[k] <- 0.5 * sum(abs(x[, k]/varTotal[k] - tx/Total))
    return(result)
}

#' A function to compute Spatial Gini's segregation index 
#'
#' @usage Gini(x)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A vector with Gini segregation index 
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' Gini(x)
#' @references Gini C. (1921) \emph{Measurement of inequality of income}. 
#' Economic Journal 31, pp. 22-43
#' @description A spatial version of Gini index. It can be derived 
#' from the Lorenz curve, and varies between 0 (minimum segregation) 
#' and 1 (maximum segregation).
#' @seealso Other evenness intragroup  indices: \code{\link{Duncan}}, 
#' \code{\link{Morill}}, \code{\link{Wong}}, \code{\link{Gorard}}
#' @seealso Intergroup dissimilarity index: \code{\link{DI}}
#' @export


Gini <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    t <- rowSums(x)
    p <- x/t
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    for (k in 1:ncol(x)) {
        for (i in 1:nrow(x)) pij[k, , i] <- p[i, k]
        for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - p[i, k])
        matprovi <- (t %*% t(t)) * pij[k, , ]
        result[k] <- sum(matprovi)/(2 * Total * Total * pTotal[k] * (1 - pTotal[k]))
    }
    return(result)
}


#' A function to compute Duncan dissimilarity segregation index
#'
#' @usage DI(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A matrix with Duncan dissimilarity index
#' @references Duncan O. D. and Duncan B. (1955) \emph{A Methodological 
#' Analysis of Segregation Indexes}. American Sociological Review 41, 
#' pp. 210-217
#' @description Duncan's dissimilarity index is the most common segregation 
#' index used to measure the segregation between different groups. For each 
#' pair of groups, the index varies between 0 (similar distribution of groups) 
#' to 1 (perfect spatial separation)
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' DI(x) 
#' @seealso Evenness intragroup  indices: \code{\link{Duncan}}, 
#' \code{\link{Wong}}, \code{\link{Morill}}, \code{\link{Gini}}, \code{\link{Gorard}}
#' @export


DI <- function(x) {
    x <- as.matrix(x)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    varTotal <- colSums(x)
    for (k1 in 1:ncol(x)) 
      for (k2 in 1:ncol(x)) 
        result[k1, k2] <- 0.5 * sum(abs(x[, k1]/varTotal[k1] - x[, k2]/varTotal[k2]))
    return(result)
}

#################### 

# EXPOSITION INDEXES

#################### 

#' A function to compute Bell's isolation index (xPx)
#'
#' @usage xPx(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A vector with isolation index 
#' @references Bell W. (1954) \emph{A probability model for the 
#' measurement of ecological segregation}. Social Forces 32(4), 
#' pp. 357-364
#' @references Massey D. S. and Denton N. A. (1988) \emph{
#' The dimensions of residential segregation}.  Social Forces 67(2),  
#' pp. 281-315.
#' @description The isolation index is a exposure index that measures 
#' the probability that a member of a group shares the same spatial 
#' unit with another member of its group.
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' xPx(x) 
#' @seealso Adjusted isolation index: \code{\link{Eta2}}
#' @seealso Interaction index: \code{\link{xPy}}
#' @export

xPx <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    t <- rowSums(x)
    for (k in 1:ncol(x)) 
      result[k] <- sum(x[, k]/varTotal[k] * x[, k]/t)
    return(result)
}


#' A function to compute adjusted isolation index (Eta2)
#'
#' @usage Eta2(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A vector with adjusted isolation index 
#' @references Bell W. (1954) \emph{A probability model for the 
#' measurement of ecological segregation}. Social Forces 32(4), 
#' pp. 357-364
#' @references Massey D. S. and Denton N. A. (1988) \emph{The 
#' dimensions of residential segregation}. Social Forces 67(2),  
#' pp. 281-315.
#' @description Adjusted isolation index (also known as Eta2 or 
#' the Correlation Ratio index) adapts isolation index \code{
#' \link{xPx}} by taking into account the differences in area 
#' population proportions.
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' Eta2(x) 
#' @seealso Bell's isolation index: \code{\link{xPx}}
#' @seealso Interaction index: \code{\link{xPy}}
#' @export

Eta2 <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    t <- rowSums(x)
    for (k in 1:ncol(x)) 
      result[k] <- sum(x[, k]/varTotal[k] * x[, k]/t)
    for (k in 1:ncol(x)) 
      result[k] <- (result[k] - pTotal[k])/(1 - pTotal[k])
    return(result)
}

#' A function to compute interaction index (xPy)
#'
#' @usage xPy(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A matrix with interaction index 
#' @references Bell W. (1954) \emph{A probability model for the 
#' measurement of ecological segregation}. Social Forces 32(4),
#'  pp. 357-364
#' @references Massey D. S. and Denton N. A. (1988) \emph{The 
#' dimensions of residential segregation}. Social Forces 67(2),  
#' pp. 281-315.
#' @description The interaction index is an exposure index that 
#' measures the probability that a member of a group shares the 
#' same spatial unit with a member of another group.
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' xPy(x) 
#' @seealso Intragroup exposure indices: \code{\link{xPx}},  
#' \code{\link{Eta2}}
#' @export


xPy <- function(x) {
    x <- as.matrix(x)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    varTotal <- colSums(x)
    t <- rowSums(x)
    for (k1 in 1:ncol(x)) 
      for (k2 in 1:ncol(x)) 
        result[k1, k2] <- sum(x[, k1]/varTotal[k1] * x[, k2]/t)
    return(result)
}


####################### 

# CONCENTATION INDEXES

####################### 


#' A function to compute Delta index
#'
#' @usage Delta(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param a - a vector with areas of spatial units
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A vector with Delta index 
#' @references Duncan O. D., Cuzzoert  and Duncan B. (1961) 
#' \emph{Problems in analyzing areal data}. Statistical geography, 
#' Glencoe, Illinois: The free press of Glencoe
#' @description Delta index compares the relative density of a 
#' spatial unit with the proportion of a group living in the 
#' same unit. The function can be used in two ways: by providing 
#' a vector with spatial units area, or a geographic source (spatial 
#' object or shape file) wich will be used to compute the area vector
#'  within the function
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' ar <- area(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' Delta(x, a = ar) 
#' 
#' Delta(x, spatobj = GreHSize)
#' 
#' Delta(x, folder = foldername, shape = shapename) 
#' @seealso Absolute Concentration Index: \code{\link{ACO}}
#' @seealso Relative Concentration Index: \code{\link{RCO}}
#' @export


Delta <- function(x, a = NULL, spatobj = NULL, folder = NULL, 
                  shape = NULL) {
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    areaTotal <- sum(a)
    for (k in 1:ncol(x)) 
      result[k] <- 0.5 * sum(abs(x[, k]/varTotal[k] - a/areaTotal))
    return(result)
}



#' A function to compute Absolute Concentration index (ACO)
#'
#' @usage ACO(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param a - a vector with areas of spatial units
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A vector with Absolute Concentration index 
#' @references Massey D. S. and Denton N. A. (1988) \emph{
#' The dimensions of residential segregation}. 
#' Social Forces 67(2),  pp. 281-315.
#' @description Concentration refers to the physical space occupied 
#' by a group. The less of the area a group occupies, the more 
#' concentrated it is. The function can be used in two ways:
#'  by providing a vector with spatial units' area, or a geographic 
#'  source (spatial object or shape file) wich will be used to 
#'  compute the area vector within the function.
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' ar <- area(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' ACO(x, a = ar) 
#' 
#' ACO(x, spatobj = GreHSize)
#' 
#' ACO(x, folder = foldername, shape = shapename) 
#' 
#' @seealso Delta Index: \code{\link{Delta}}
#' @seealso Relative Concentration Index: \code{\link{RCO}}
#' @export


ACO <- function(x, a = NULL, spatobj = NULL, folder = NULL, 
                shape = NULL) {
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    varTotal <- colSums(x)
    xprovi <- as.data.frame(cbind(x, a))
    xprovi <- xprovi[order(xprovi$a), ]
    xprovi$Total <- rowSums(xprovi) - xprovi$a
    areaTotal <- sum(a)
    result <- vector(length = ncol(x))
    n1 <- vector(length = ncol(x))
    n2 <- vector(length = ncol(x))
    T1 <- vector(length = ncol(x))
    T2 <- vector(length = ncol(x))
    t <- rowSums(xprovi)
    
    for (k in 1:ncol(x)) {
        T1[k] <- 0
        i <- 0
        while (T1[k] < varTotal[k]) {
            i <- i + 1
            T1[k] <- T1[k] + xprovi$Total[i]
        }
        n1[k] <- i
        T2[k] <- 0
        i <- nrow(xprovi) + 1
        while (T2[k] < varTotal[k]) {
            i <- i - 1
            T2[k] <- T2[k] + xprovi$Total[i]
        }
        n2[k] <- i
    }
    
    for (k in 1:ncol(x)) {
        vartemp1 <- sum(xprovi[, k] * xprovi$a/varTotal[k])
        vartemp2 <- sum(xprovi$Total[1:n1[k]] * xprovi$a[1:n1[k]]/T1[k])
        vartemp3 <- sum(xprovi$Total[n2[k]:nrow(xprovi)] * xprovi$a[n2[k]:nrow(xprovi)]/T2[k])
        result[k] <- 1 - (vartemp1 - vartemp2)/(vartemp3 - vartemp2)
    }
    
    return(result)
}


#' A function to compute Relative Concentration index (RCO)
#'
#' @usage RCO(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param a - a vector with areas of spatial units
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A matrix with relative concentration index
#' @references Massey D. S. and Denton N. A. (1988) \emph{
#' The dimensions of residential segregation}. 
#' Social Forces 67(2),  pp. 281-315.
#' @description Relative concentration index measures the share 
#' of urban space occupied by a group compared to another group. 
#' The function can be used in two ways: by providing a vector with 
#' spatial units' area or a geographic source (spatial object or 
#' shape file) wich will be used to compute the area vector within 
#' the function
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' ar <- area(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' RCO(x, a = ar) 
#' 
#' RCO(x, spatobj = GreHSize)
#' 
#' RCO(x, folder = foldername, shape = shapename) 
#' 
#' @seealso Intragroup concentration indices: \code{\link{Delta}},  \code{\link{ACO}}
#' @export


RCO <- function(x, a = NULL, spatobj = NULL, folder = NULL, 
                shape = NULL) {
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    varTotal <- colSums(x)
    xprovi <- as.data.frame(cbind(x, a))
    xprovi <- xprovi[order(xprovi$a), ]
    xprovi$Total <- rowSums(xprovi) - xprovi$a
    areaTotal <- sum(a)
    n1 <- vector(length = ncol(x))
    n2 <- vector(length = ncol(x))
    T1 <- vector(length = ncol(x))
    T2 <- vector(length = ncol(x))
    t <- rowSums(xprovi)
    
    for (k in 1:ncol(x)) {
        T1[k] <- 0
        i <- 0
        while (T1[k] < varTotal[k]) {
            i <- i + 1
            T1[k] <- T1[k] + xprovi$Total[i]
        }
        n1[k] <- i
        T2[k] <- 0
        i <- nrow(xprovi) + 1
        while (T2[k] < varTotal[k]) {
            i <- i - 1
            T2[k] <- T2[k] + xprovi$Total[i]
        }
        n2[k] <- i
    }
    for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) {
        vartemp1 <- sum(xprovi[, k1] * xprovi$a/varTotal[k1])
        vartemp1b <- sum(xprovi[, k2] * xprovi$a/varTotal[k2])
        vartemp2 <- sum(xprovi$Total[1:n1[k1]] * xprovi$a[1:n1[k1]]/T1[k1])
        vartemp3 <- sum(xprovi$Total[n2[k1]:nrow(xprovi)] * xprovi$a[n2[k1]:nrow(xprovi)]/T2[k1])
        result[k1, k2] <- ((vartemp1/vartemp1b) - 1)/((vartemp2/vartemp3) - 1)
    }
    return(result)
}

####################### 

# CLUSTERING INDEXES

####################### 


#' A function to compute Absolute Clustering Index (ACL)
#'
#' @usage ACL(x, c = NULL, spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param c - a binary contiguity (adjacency) symetric matrix where each 
#' element \emph{Cij} equals 1 if \emph{i}-th and \emph{j}-th spatial units 
#' are adjacent, 0 otherwise. 
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A vector with Absolute Clustering index 
#' @references Massey D. S. and Denton N. A. (1988) \emph{
#' The dimensions of residential segregation}. 
#' Social Forces 67(2),  pp. 281-315.
#' @description The more contiguous spatial units a group occupies 
#' (forming an enclave within the zone) the more clustered and therefore 
#' segregated it is. The function can be used in two ways: 
#'  by providing a contiguity matrix or a geographic source (spatial 
#'  object or shape file) wich will be used to compute the matrix within 
#'  the function
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' contiguity <- contig(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' ACL(x, c=contiguity) 
#' 
#' ACL(x, spatobj = GreHSize)
#' 
#' ACL(x, folder = foldername, shape = shapename) 
#'
#' @seealso Mean proximity between members of a group: \code{\link{Pxx}}
#' @seealso Relative Clustering Index: \code{\link{RCL}}
#' @export


ACL <- function(x, c = NULL, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(c)) 
        c <- contig(spatobj = spatobj, folder = folder, shape = shape)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    for (k in 1:ncol(x)) {
        vartemp1 <- sum((c %*% x[, k]) * x[, k]/varTotal[k])
        t <- as.vector(rowSums(x))
        vartemp2 <- (sum(c) * varTotal[k])/(nrow(x)^2)
        vartemp3 <- sum((c %*% t) * x[, k]/varTotal[k])
        result[k] <- (vartemp1 - vartemp2)/(vartemp3 - vartemp2)
    }
    return(result)
}


#' A function to compute the mean proximity between members of a group (Pxx)
#'
#' @usage Pxx(x, d = NULL, fdist = 'l', spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param d - a matrix with the distances between spatial units centroids
#' @param fdist - the method used for distance calculations: 
#' 'l' for linear (by default) and 'e' for exponential function.
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A vector with Pxx values
#' @references White M. J. (1988) \emph{The Measurement of Spatial 
#' Segregation}. American Journal of Sociology, 88, p. 1008-1019
#' @description The index computes the mean distance between the members
#'  of a group. The function can be used in two ways: by providing a 
#'  distance matrix or a geographic source (spatial object or shape file) 
#'  wich will be used to compute the matrix within the function
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' dist <- distance(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' Pxx(x, d = dist) 
#' 
#' Pxx(x, spatobj = GreHSize) 
#' 
#' Pxx(x, folder = foldername, shape = shapename) 
#' 
#' @seealso Absolute Clustering Index: \code{\link{ACL}}
#' @seealso Other proximity measures: \code{\link{Poo}}, 
#' \code{\link{Pxy}}, \code{\link{SP}}
#' @export


Pxx <- function(x, d = NULL, fdist = "l", spatobj = NULL, 
                folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape)
    varTotal <- colSums(x)
    d2 <- d
    if (fdist == "e") 
        d2 <- exp(-d)
    result <- vector(length = ncol(x))
    for (k in 1:ncol(x)) 
      result[k] <- sum((d2 %*% x[, k]) * x[, k]/(varTotal[k])^2)
    return(result)
}

#' A function to compute the mean proximity between 
#' persons without regard to group (Poo)
#'
#' @usage Poo(x, d = NULL, fdist = 'l', spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param d - a matrix with the distances between spatial 
#' units centroids
#' @param fdist - the method used for distance calculations: 
#' 'l' for linear (by default) and 'e' for exponential function.
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A matrix with Poo values
#' @references White M. J. (1988) \emph{The Measurement of Spatial 
#' Segregation}. American Journal of Sociology, 88, p. 1008-1019
#' @description The index computes the mean distance between the 
#' persons in the area without regard to group. The function can 
#' be used in two ways: by providing a distance matrix or a 
#' geographic source (spatial object or shape file) wich will be 
#' used to compute the matrix within the function
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' dist<- distance(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' Poo(x, d = dist) 
#' 
#' Poo(x, spatobj = GreHSize) 
#' 
#' Poo(x, folder = foldername, shape = shapename) 
#'
#' @seealso Absolute Clustering Index: \code{\link{ACL}}
#' @seealso Other proximity measures: \code{\link{Pxx}} ,
#' \code{\link{Pxy}}, \code{\link{SP}}
#' @export


Poo <- function(x, d = NULL, fdist = "l", spatobj = NULL, 
                folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape)
    varTotal <- colSums(x)
    d2 <- d
    if (fdist == "e") 
        d2 <- exp(-d)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    for (k1 in 1:ncol(x)) 
      for (k2 in 1:ncol(x)) 
        result[k1, k2] <- sum((d2 %*% (x[, k1] + x[, k2])) * (x[, k1] + 
        x[, k2])/(varTotal[k1] + varTotal[k2])^2)
    return(result)
}

#' A function to compute the mean proximity between 
#' persons of different groups (Pxy)
#'
#' @usage Pxy(x, d = NULL, fdist = 'l', spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param d - a matrix with the distances between spatial 
#' units centroids
#' @param fdist - the method used for distance calculations: 
#' 'l' for linear (by default) and 'e' for exponential function.
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A matrix with Pxx values
#' @references White M. J. (1988) \emph{The Measurement of 
#' Spatial Segregation}. 
#' American Journal of Sociology, 88, p. 1008-1019
#' @description The Pxy index computes the mean distance between 
#' the persons of different groups.  The function can be used in 
#' two ways: by providing a distance matrix, or a geographic source 
#' (spatial object or shape file) wich will be used to compute the 
#' matrix within the function.
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' dist<- distance(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' Pxy(x, d = dist) 
#' 
#' Pxy(x, spatobj = GreHSize) 
#' 
#' Pxy(x, folder = foldername, shape = shapename) 
#'
#' @seealso Absolute Clustering Index: \code{\link{ACL}}
#' @seealso Other proximity measures: \code{\link{Pxx}} , 
#' \code{\link{Poo}}, \code{\link{SP}}
#' @export


Pxy <- function(x, d = NULL, fdist = "l", spatobj = NULL, 
                folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape)
    varTotal <- colSums(x)
    d2 <- d
    if (fdist == "e") 
        d2 <- exp(-d)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    for (k1 in 1:ncol(x)) 
      for (k2 in 1:ncol(x)) 
        result[k1, k2] <- sum((d2 %*% (x[, k2])) * (x[, k1]))/varTotal[k1]/varTotal[k2]
    return(result)
}

#' A function to compute the proximity index (SP)
#'
#' @usage SP(x, d = NULL, fdist = 'l', spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param d - a matrix with the distances between spatial 
#' units centroids
#' @param fdist - the method used for distance calculations: 
#' 'l' for linear (by default) and 'e' for exponential function.
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A matrix with SP values
#' @references White M. J. (1988) \emph{The Measurement of Spatial 
#' Segregation}. American Journal of Sociology, 88, p. 1008-1019
#' @description The index compares the clustering level of a group 
#' to another. The function can be used in two ways: by providing 
#' a distance matrix or a geographic source (spatial object or shape 
#' file) wich will be used to compute the matrix within the function
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' dist <- distance(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' SP(x, d = dist) 
#' 
#' SP(x, spatobj = GreHSize) 
#' 
#' SP(x, folder = foldername, shape = shapename) 
#'
#' @seealso Absolute Clustering Index: \code{\link{ACL}}
#' @seealso Other proximity measures: \code{\link{Pxx}}, \code{\link{Pxy}}, 
#' \code{\link{Poo}}
#' @export


SP <- function(x, d = NULL, fdist = "l", spatobj = NULL, folder = NULL, 
               shape = NULL) {
    x <- as.matrix(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape)
    varTotal <- colSums(x)
    d2 <- d
    if (fdist == "e") 
        d2 <- exp(-d)
    PxxExp <- vector(length = ncol(x))
    for (k in 1:ncol(x)) 
      PxxExp[k] <- sum((d2 %*% x[, k]) * x[, k]/(varTotal[k])^2)
    PooExp <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    for (k1 in 1:ncol(x)) 
      for (k2 in 1:ncol(x)) 
        PooExp[k1, k2] <- sum((d2 %*% (x[, k1] + x[, k2])) * (x[, k1] + x[, k2])/
                                (varTotal[k1] + varTotal[k2])^2)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    for (k1 in 1:ncol(x)) 
      for (k2 in 1:ncol(x)) 
        result[k1, k2] <- (varTotal[k1] * PxxExp[k1] + varTotal[k2] * PxxExp[k2])/
        ((varTotal[k1] + varTotal[k2]) * PooExp[k1, k2])
    return(result)
}



#' A function to compute the relative clustering index (RCL)
#'
#' @usage RCL(x, d = NULL, fdist = 'l', spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param d - a matrix with the distances between spatial units 
#' centroids
#' @param fdist - the method used for distance calculations: 
#' 'l' for linear (by default) and 'e' for exponential function.
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A matrix with RCL values
#' @references Massey D. S. and Denton N. A. (1988) \emph{The dimensions 
#' of residential segregation}. Social Forces 67(2),  pp. 281-315.
#' @description The index compares the mean proximity of a group and 
#' its mean promiximy to another group. The function can be used in 
#' two ways: by providing a distance matrix or a geographic source 
#' (spatial object or shape file) wich will be used to compute the
#'  matrix within the function
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' dist<- distance(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' RCL(x, d = dist) 
#' 
#' RCL(x, spatobj = GreHSize) 
#' 
#' RCL(x, folder = foldername, shape = shapename) 
#' @seealso Absolute Clustering Index: \code{\link{ACL}}
#' @seealso Other intergroup clustering measures:  \code{\link{Pxy}}, 
#' \code{\link{SP}}
#' @export


RCL <- function(x, d = NULL, fdist = "l", spatobj = NULL, 
                folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape)
    varTotal <- colSums(x)
    d2 <- exp(-d)
    PxxExp <- vector(length = ncol(x))
    for (k in 1:ncol(x)) 
      PxxExp[k] <- sum((d2 %*% x[, k]) * x[, k]/(varTotal[k])^2)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    for (k1 in 1:ncol(x)) 
      for (k2 in 1:ncol(x)) 
        result[k1, k2] <- PxxExp[k1]/PxxExp[k2] - 1
    return(result)
}

######################## 

# CENTRALISATION INDEXES

######################## 

#' A function to compute Absolute Centralisation Index (ACE)
#'
#' @usage ACE(x, a = NULL, dc = NULL, center = 1, spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param a - a vector with areas of spatial units
#' @param dc - a vector with the distances between spatial units 
#' centroids and the center
#' @param center - a value giving the number of the spatial unit 
#' that represents the zone's center
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A vector with Absolute Centralisation index
#' @references Duncan O. D. and Duncan B. (1955) \emph{A 
#' Methodological Analysis of Segregation Indexes}. 
#' American Sociological Review 41, pp. 210-217
#' @description Absolute centralisation index measures the proportion 
#' of a group that should change localisation to obtain a uniform 
#' density around the area's center. The function can be used in 
#' two ways: by providing  the area vector and the distances to the 
#' center vector or a geographic source (spatial object or shape file) 
#' wich will be used to compute the vectors within the function. 
#' The center parameter is necessary to specify the number of the spatial 
#' unit representing the center
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' ar<-area(GreHSize)
#' distc<- distcenter(GreHSize, center = 19)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' ACE(x, a = ar, dc=distc) 
#' 
#' ACE(x, spatobj = GreHSize, center = 19) 
#' 
#' ACE(x, folder = foldername, shape = shapename, center = 19) 
#'
#' @seealso Relative Centralisation Index: \code{\link{RCE}}
#' @export


ACE <- function(x, a = NULL, dc = NULL, center = 1, 
                spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    if (is.null(dc)) 
        dc <- distcenter(spatobj = spatobj, folder = folder, 
                         shape = shape, center)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    xprovi <- cbind(x, a, dc)
    xprovi <- xprovi[order(xprovi[, ncol(xprovi)]), ]
    xprovi <- as.data.frame(xprovi)
    for (k in 1:ncol(x)) {
        XI1 <- cumsum(xprovi[, k])[1:(nrow(xprovi) - 1)]/varTotal[k]
        AI <- cumsum(xprovi$a)[2:nrow(xprovi)]/sum(xprovi$a)
        XI <- cumsum(xprovi[, k])[2:nrow(xprovi)]/varTotal[k]
        AI1 <- cumsum(xprovi$a)[1:(nrow(xprovi) - 1)]/sum(xprovi$a)
        result[k] <- XI1 %*% AI - XI %*% AI1
    }
    return(result)
}



#' A function to compute Relatice Centralisation Index (RCE)
#'
#' @usage RCE(x, dc = NULL, center = 1, spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param dc - a vector with distances from the spatial units cetroids 
#' to the center
#' @param center - a value giving the number of the spatial unit 
#' that represents the zones's center
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A matrix with Relative Centralisation index 
#' @references Duncan O. D. and Duncan B. (1955) \emph{A 
#' Methodological Analysis of Segregation Indexes}. 
#' American Sociological Review 41, pp. 210-217
#' @description Relative Centralisation Index measures the proportion 
#' of a group that should change localisation to obtain the same level 
#' of centralisation that another group. The function can be used 
#' in two ways: by providing a vector with the distances to the center 
#' or a geographic source (spatial object or shape file) wich will be 
#' used to compute the vector within the function. The center parameter 
#' is necessary to specify the number of the spatial unit representing 
#' the center.
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' distc<- distcenter(GreHSize, center = 19)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' RCE(x, dc=distc) 
#' 
#' RCE(x, spatobj = GreHSize, center = 19) 
#' 
#' RCE(x, folder = foldername, shape = shapename, center = 19) 
#'
#' @seealso Absolute Centralisation Index: \code{\link{ACE}}
#' @export


RCE <- function(x, dc = NULL, center = 1, spatobj = NULL, 
                folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(dc)) 
        dc <- distcenter(spatobj = spatobj, folder = folder, 
                         shape = shape, center)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    varTotal <- colSums(x)
    xprovi <- cbind(x, dc)
    xprovi <- xprovi[order(xprovi[, ncol(xprovi)]), ]
    xprovi <- as.data.frame(xprovi)
    for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) {
        XI1 <- cumsum(xprovi[, k1])[1:(nrow(xprovi) - 1)]/varTotal[k1]
        XI <- cumsum(xprovi[, k1])[2:nrow(xprovi)]/varTotal[k1]
        YI1 <- cumsum(xprovi[, k2])[1:(nrow(xprovi) - 1)]/varTotal[k2]
        YI <- cumsum(xprovi[, k2])[2:nrow(xprovi)]/varTotal[k2]
        result[k1, k2] <- XI1 %*% YI - XI %*% YI1
    }
    return(result)
}


######################## 

# MULTIGROUP INDEXES

######################## 


#' A function to compute Shannon-Wiener diversity index
#'
#' @usage HShannon(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return Shannon-Wiener diversity index 
#' @references Shannon C. E. (1948) \emph{A mathematical theory 
#' of communication}. Bell System Technical Journal (27) 
#' @description Shannon-Wiener diversity index is based on entropy 
#' notion and measures the population heterogeneity. 
#' @examples  x <- slot(GreHSize, 'data')[ ,3:5]
#' HShannon(x) 
#' @seealso  Other multigroup eveness indices: 
#' \code{\link{JPielou}}, \code{\link{ISimpson}}, 
#' \code{\link{GiniMulti}}, \code{\link{DMulti}}
#' @seealso Other multigroup indices: \code{\link{PIsol}}, 
#' \code{\link{RelDivers}}
#' @export


HShannon <- function(x) {
    x <- as.matrix(x)
    pTotal <- colSums(x)/sum(x)
    result <- -sum(pTotal * log(pTotal))
    return(result)
}


#' A function to compute Pielou normalised diversity index
#'
#' @usage JPielou(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return Normalised diversity index 
#' @references Pielou E. C. (1966) \emph{Shannon's formula as a 
#' measure of species diversity: its use and misuse}. The 
#' American Naturalist 100, 463-465 
#' @description Pielou index normalises Shannon-Wiener diversity 
#' index (\code{\link{HShannon}}).
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' JPielou(x) 
#' @seealso  Other multigroup eveness indices: 
#' \code{\link{HShannon}}, \code{\link{ISimpson}}, \code{\link{GiniMulti}}, 
#' \code{\link{DMulti}}
#' @seealso Other multigroup indices: \code{\link{PIsol}}, 
#' \code{\link{RelDivers}}
#' @export


JPielou <- function(x) {
    x <- as.matrix(x)
    pTotal <- colSums(x)/sum(x)
    result <- -sum(pTotal * log(pTotal))
    result <- result/log(ncol(x))
    return(result)
}


#' A function to compute Simpson's interaction index
#'
#' @usage ISimpson(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return Simpson's interaction index 
#' @references Simpson E. H. (1949) \emph{Measurement of diversity}. 
#' Nature 163:688 
#' @description Simpson's interaction index measures the probability 
#' that individuals randomly selected are not the same group.
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' ISimpson(x) 
#' @seealso  Other multigroup eveness indices: 
#' \code{\link{HShannon}}, \code{\link{JPielou}}, \code{\link{GiniMulti}}, 
#' \code{\link{DMulti}}
#' @seealso Other multigroup indices: \code{\link{PIsol}}, 
#' \code{\link{RelDivers}}
#' @export


ISimpson <- function(x) {
    x <- as.matrix(x)
    pTotal <- colSums(x)/sum(x)
    result <- sum(pTotal * (1 - pTotal))
    return(result)
}



#' A function to compute multigroup Gini index
#'
#' @usage GiniMulti(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @references Reardon S. F. (1998) \emph{Measures of racial 
#' diversity and segregation in multigroup and hierarchical 
#' structured Populations}. Annual meeting of the Eastern 
#' Sociological Society, Philadelphia 
#' @description Multigroup version of \code{\link{Gini}} index 
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' GiniMulti(x) 
#' @seealso  Other multigroup eveness indices: 
#' \code{\link{HShannon}}, \code{\link{JPielou}}, 
#' \code{\link{ISimpson}}, \code{\link{DMulti}}
#' @seealso Other multigroup indices: \code{\link{PIsol}}, 
#' \code{\link{RelDivers}}
#' @export


GiniMulti <- function(x) {
    x <- as.matrix(x)
    pTotal <- colSums(x)/sum(x)
    II <- sum(pTotal * (1 - pTotal))
    vartemp <- vector(length = ncol(x))
    t <- rowSums(x)
    Total <- sum(x)
    p <- x/t
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    for (k in 1:ncol(x)) {
        for (i in 1:nrow(x)) pij[k, , i] <- p[i, k]
        for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - p[i, k])
        vartemp[k] <- sum(t * (t %*% pij[k, , ]))
    }
    result <- sum(vartemp)/(2 * Total * Total * II)
    return(result)
}





#' A function to compute multigroup dissimilarity index
#'
#' @usage DMulti(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return Multigroup dissimilarity index 
#' @references Sakoda J. N. (1981) \emph{A generalized Index of 
#' dissimilarity}. Demography,18, 245-250 
#' @description Multigroup version of dissimilarity index 
#' (\code{\link{DI}})
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' DMulti(x) 
#' @seealso  Other multigroup eveness indices: 
#' \code{\link{HShannon}}, \code{\link{JPielou}}, \code{\link{ISimpson}}, 
#' \code{\link{GiniMulti}}
#' @seealso Other multigroup indices: \code{\link{PIsol}}, 
#' \code{\link{RelDivers}}
#' @export


DMulti <- function(x) {
    x <- as.matrix(x)
    result <- 0
    pTotal <- colSums(x)/sum(x)
    II <- sum(pTotal * (1 - pTotal))
    t <- rowSums(x)
    Total <- sum(x)
    p <- x/t
    for (k in 1:ncol(x)) 
      result <- result + t %*% abs(p[, k] - pTotal[k])
    result <- result/(2 * Total * II)
    return(result)
}




#' A function to compute multigroup normalised isolation index (PIsol)
#'
#' @usage PIsol(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return Multigroup normalised isolation index 
#' @references Massey D. S. and Denton N. A. (1988) \emph{The 
#' dimensions of residential segregation}. 
#' Social Forces 67(2),  pp. 281-315.
#' @description Multigroup version of isolation index (\code{\link{xPx}})
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' PIsol(x) 
#' @seealso Other multigroup exposition indices:  \code{\link{RelDivers}}
#' @seealso  Other multigroup indices: 
#' \code{\link{HShannon}}, \code{\link{JPielou}}, \code{\link{ISimpson}}, 
#' \code{\link{GiniMulti}}, \code{\link{DMulti}}
#' @export


PIsol <- function(x) {
    x <- as.matrix(x)
    result <- 0
    pTotal <- colSums(x)/sum(x)
    t <- rowSums(x)
    Total <- sum(x)
    p <- x/t
    for (k in 1:ncol(x)) 
      result <- result + t %*% ((p[, k] - pTotal[k]) * 
                                  (p[, k] - pTotal[k])/(1 - p[, k]))
    result <- result/Total
    return(result)
}


#' A function to compute multigroup relative diversity index
#'
#' @usage RelDivers(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return Multigroup relative diversity index 
#' @references Carlson S. M. (1992) \emph{Trends in race/sex 
#' occupational inequality:  conceptual and measurement issues}. 
#' Social Problems, 39, p. 269-290
#' @description Multigroup index based on Simpson's interaction 
#' index \code{\link{ISimpson}}
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' RelDivers(x) 
#' @seealso Other multigroup exposition indices:  \code{\link{PIsol}}
#' @seealso Other multigroup indices: \code{\link{HShannon}}, 
#' \code{\link{JPielou}}, \code{\link{ISimpson}}, \code{\link{GiniMulti}}, 
#' \code{\link{DMulti}}
#' @export


RelDivers <- function(x) {
    x <- as.matrix(x)
    result <- 0
    pTotal <- colSums(x)/sum(x)
    II <- sum(pTotal * (1 - pTotal))
    t <- rowSums(x)
    Total <- sum(x)
    p <- x/t
    for (k in 1:ncol(x)) 
      for (k in 1:ncol(x)) 
        result <- result + t %*% ((p[, k] - pTotal[k]) * (p[, k] - pTotal[k]))
    result <- result/(Total * II)
    return(result)
}


######################## 

# LOCAL INDEXES

######################## 

#' A function to compute location quotients (LQs)
#'
#' @usage LQ(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A matrix with LQs
#' @references Isard W. (1960) \emph{Methods of regional analysis: 
#' an introduction to regional science}. The MIT Press, Cambridge
#' @description Location quotients compare the relative part of 
#' group in each spatial unit to the relative part of that same 
#' group in the area.
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' LQ(x) 
#' @seealso Other local indices:  local entropy index \code{\link{HLoc}}
#' @export

LQ <- function(x) {
    x <- as.matrix(x)
    result <- x
    Total <- sum(x)
    t <- rowSums(x)
    for (i in 1:ncol(x)) 
      result[, i] <- (x[, i]/t)/(sum(x[, i])/sum(t))
    return(result)
}


#' A function to compute local diversity (entropy) index
#'
#' @usage HLoc(x) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @return A vector with local diversity (entropy) index
#' @references Theil H. (1972) \emph{Statistical Decomposition Analysis}. 
#' North-Holland, Amsterdam
#' @description Local adaptation of Pielou's normalised diversity 
#' index \code{\link{JPielou}}
#' @examples x <- slot(GreHSize, 'data')[ ,3:5]
#' HLoc(x) 
#' @seealso Other local indices:  location quotients \code{\link{LQ}}
#' @export


HLoc <- function(x) {
    x <- as.matrix(x)
    Total <- sum(x)
    t <- rowSums(x)
    result <- x/t
    result <- cbind(result, 0)
    for (i in 1:nrow(result)) {
        n <- 0
        for (j in 1:(ncol(result) - 1)) 
          if (result[i, j] > 0) {
            result[i, ncol(result)] <- result[i, ncol(result)] + result[i, j] * log(result[i, j])
            n <- n + 1
            }
        result[i, ncol(result)] <- -result[i, ncol(result)]/log(n)
    }
    result <- result[, ncol(result)]
    return(result)
}


######################## 

# XSEG FUNCTION

######################## 



#' A function to compute all segregation indices developped in OasisR
#'
#' @usage xseg(x, center = 1, fdist = 'l', spatobj = NULL, folder = NULL, shape = NULL ) 
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param center - a value giving the number of the spatial unit 
#'  that represents the area's center
#' @param fdist - the method used for distance calculations: 
#'  'l' for linear (by default) and 'e' for exponential function.
#' @param folder - a character vector with the folder (directory) 
#' where the shapefile is
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A list with segregation indices 
#' @references Tivadar M., Schaeffer Y, Torre A. and Bray F. (2014) \emph{OASIS - 
#' un Outil d'Analyse de la Segregation et des Inegalites Spatiales}. 
#' Cybergeo : European Journal of Geography, GeOpenMod, document 699 
#' @description A general function that computes all segregation indexes available in OasisR package
#' @examples x <- slot(AnnHAge, 'data')[ ,3:5]
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'AnnHAge'
#' 
#' xseg(x, center = 1, spatobj = AnnHAge)
#' 
#' xseg(x, center = 1, folder = foldername, shape = shapename)
#' 
#' @seealso \code{\link{xgeo}} function which computes all geografical information necessary for segregation indices
#' @seealso \code{\link{MCTest}} function which computes Monte Carlos simulations to test segregation indexes
#' @export

xseg <- function(x, center = 1, fdist = "l", spatobj = NULL, folder = NULL, shape = NULL) {
    geoinfo <- xgeo(spatobj = spatobj, folder = folder, shape = shape, center = center)
    vDuncan <- Duncan(x)
    vMorill <- Morill(x, c = geoinfo$contiguity)
    vWong <- Wong(x, b = geoinfo$boundaries, p = geoinfo$perimeter, a = geoinfo$area)
    vGorard <- Gorard(x)
    vGini <- Gini(x)
    vDI <- DI(x)
    vxPx <- xPx(x)
    vEta2 <- Eta2(x)
    vxPy <- xPy(x)
    vDelta <- Delta(x, a = geoinfo$area)
    vACO <- ACO(x, a = geoinfo$area)
    vRCO <- RCO(x, a = geoinfo$area)
    vACL <- ACL(x, c = geoinfo$contiguity)
    vPxx <- Pxx(x, d = geoinfo$distance, fdist = fdist)
    vPoo <- Poo(x, d = geoinfo$distance, fdist = fdist)
    vPxy <- Pxy(x, d = geoinfo$distance, fdist = fdist)
    vSP <- SP(x, d = geoinfo$distance, fdist = fdist)
    vRCL <- RCL(x, d = geoinfo$distance, fdist = fdist)
    vACE <- ACE(x, a = geoinfo$area, dc = geoinfo$distcenter)
    vRCE <- RCE(x, dc = geoinfo$distcenter)
    vHShannon <- HShannon(x)
    vJPielou <- JPielou(x)
    vISimpson <- ISimpson(x)
    vGiniMulti <- GiniMulti(x)
    vDMulti <- DMulti(x)
    vPIsol <- PIsol(x)
    vRelDivers <- RelDivers(x)
    vLQ <- LQ(x)
    vHLoc <- HLoc(x)
    resultlist <- list(Duncan = vDuncan, Morill = vMorill, Wong = vWong, Gorard = vGorard, Gini = vGini, DI = vDI, 
        xPx = vxPx, Eta2 = vEta2, xPy = vxPy, Delta = vDelta, ACO = vACO, RCO = vRCO, ACL = vACL, Pxx = vPxx, Poo = vPoo, 
        Pxy = vPxy, SP = vSP, RCL = vRCL, ACE = vACE, RCE = vRCE, HShannon = vHShannon, JPielou = vJPielou, ISimpson = vISimpson, 
        GiniMulti = vGiniMulti, DMulti = vDMulti, PIsol = vPIsol, RelDivers = vRelDivers, LQ = vLQ, HLoc = vHLoc)
    return(resultlist)
}


################################################## 

# MONTE CARLO SIMULATIONS

################################################## 


#' A function to test segregation indexes using Monte Carlo simulations
#'
#' @usage MCTest(x, fun, type = 'user', proba = NULL, c = NULL,
#'  b = NULL, p=NULL, a = NULL, d = NULL, dc = NULL, fdist = 'l', 
#'  center = 1, nsim = 99, spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or that can be coerced 
#' to that class), where each column represents the distribution 
#' of a population group, within spatial units. The number of 
#' columns should be greater than one (at least two population 
#' groups are required).
#' @param fun - a character vector with the segregation function 
#' to be tested (only intragroup or multigroup indexes)
#' @param type - a character vector with type of simulation. 
#' If type='perm', the function will use the permutation test. 
#' If type='rand', the population is realocated in the spatial 
#' units randomly. If type='area', the location probabilities 
#' are proportional to the spatial units' area. For the default 
#' value type='user', the function expects you to introduce a 
#' probability location vector.
#' @param proba - a vector with location probabilities. By 
#' default proba = NULL being calculated depending on the 
#' simulation type. The parameter is necessary when a user method 
#' is specified.
#' @param c - a binary contiguity (adjacency) symetric matrix 
#' where each element \emph{Cij} equals 1 if \emph{i}-th and 
#' \emph{j}-th spatial units are adjacent, 0 otherwise. This 
#' parameter is necessary only for some segregation functions.
#' @param b - a common boundaries matrix where each element 
#' \emph{Bij} equals  the shared boundary of \emph{i}-th and 
#' \emph{j}-th spatial units. This parameter is necessary only 
#' for some segregation functions.
#' @param p - a vector with the perimeters of spatial units
#' @param a - a vector with areas of spatial units
#' @param d - a matrix with the distances between patial units 
#' centroid
#' @param dc - a vector with distances from the spatial units 
#' cetroids to the center
#' @param fdist - a parameter with the method used for distance 
#' calculations: 'l' for linear (by default) and 'e' for 
#' exponential function.
#' @param center - a value giving the number of the spatial unit 
#' that represents the area's center
#' @param nsim - the number of simulations
#' @param folder - a character vector with the path where the shape 
#' file is located
#' @param shape - a character vector with the shape file name
#'  #' @param spatobj - a spatial object (SpatialPolygonsDataFrame)
#' @return A data frame with the simulations mean value of the index, 
#' the rank of the index in the  simulation distribution and its pseudo 
#' p.value
#' @references Tivadar M., Schaeffer Y, Torre A. and Bray F. (2014) 
#' \emph{OASIS - un Outil d'Analyse de la Segregation et des Inegalites 
#' Spatiales}.  Cybergeo : European Journal of Geography, GeOpenMod, 
#' document 699
#' @description Monte Carlo tests for intragroup or multigroup 
#' segregation indexes. 
#' @examples x <- slot(AnnHAge, 'data')[ ,3:5]
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'AnnHAge'
#' ar <- area(AnnHAge)
#' distc<- distcenter(AnnHAge, center = 19)
#' probavector<-ar/sum(ar)
#' 
#' MCTest(x, fun='Duncan', type='rand', nsim=999)
#' 
#' MCTest(x, fun='ACE', type='perm', a = ar , dc=distc)
#' 
#' MCTest(x, fun='Wong', type='area', folder = foldername, 
#' shape = shapename)
#' 
#' MCTest(x, fun='Morill', type='user', proba=probavector, 
#' spatobj = AnnHAge)
#' 
#' @seealso \code{\link{xgeo}} function which computes all geografical 
#' information necessary for segregation indices
#' @export


MCTest <- function(x, fun, type = "user", proba = NULL, c = NULL, 
                   b = NULL, p = NULL, a = NULL, d = NULL, dc = NULL, 
                   fdist = "l", center = 1, nsim = 99, spatobj = NULL, 
                   folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    assign("func", eval(parse(text = fun)))
    xdistr <- vector("list", nsim)
    
    if (type == "perm") 
        for (k in 1:nsim) {
            xdistr[[k]] <- matrix(nrow = nrow(x), ncol = ncol(x))
            neworder <- sample(c(1:nrow(x)), size = nrow(x))
            for (i in 1:nrow(x)) 
              for (j in 1:ncol(x)) 
                xdistr[[k]][i, j] <- x[neworder[i], j]
        }
    if (type == "rand") 
        proba <- rep(1/nrow(x), nrow(x))
    if (type == "area") {
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        proba <- a/sum(a)
    }
    
    if (type != "perm") {
        xprovi <- NULL
        for (k in 1:ncol(x)) 
          for (i in 1:(nrow(x))) 
            xprovi <- c(xprovi, rnorm(nsim, mean = sum(x[, k]) * proba[i], 
        sd = (sum(x[, k]) * proba[i] * (1 - proba[i]))^0.5))
        for (j in 1:length(xprovi)) 
          if (xprovi[j] < 0) 
            xprovi[j] <- 0
        dim(xprovi) <- c(nsim, nrow(x), ncol(x))
        for (i in 1:nsim) 
          for (k in 1:ncol(x)) 
            xprovi[i, , k] <- sum(x[, k]) * (xprovi[i, , k]/sum(xprovi[i, , k]))
        for (i in 1:nsim) 
          xdistr[[i]] <- xprovi[i, , ]
    }
    
    nvar <- ncol(x)
    if (fun == "HShannon" || fun == "JPielou" || fun == "ISimpson" || 
          fun == "GiniMulti" || fun == "DMulti" || fun == "PIsol" || 
          fun == "RelDivers") 
        nvar <- 1
    IndTest <- matrix(nrow = nvar, ncol = 5)
    IndTest <- as.data.frame(IndTest)
    names(IndTest) <- c("Var", fun, "Simulated", "Rank", "P.Value")
    IndTest$Var <- 1:nvar
    
    resim <- matrix(nrow = nvar, ncol = nsim)
    if (fun == "Duncan" || fun == "Gorard" || fun == "Gini" || 
          fun == "xPx" || fun == "Eta2") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]])
        IndTest[, 2] <- func(x)
    }
    if (fun == "Morill" || fun == "ACL") {
        if (is.null(c)) 
            c <- contig(spatobj = spatobj, folder = folder, shape = shape)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], c)
        IndTest[, 2] <- func(x, c)
    }
    if (fun == "Wong") {
        if (is.null(b)) 
            b <- boundaries(spatobj = spatobj, folder = folder, shape = shape)
        if (is.null(p)) 
            p <- perimeter(spatobj = spatobj, folder = folder, shape = shape)
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], b, p, a)
        IndTest[, 2] <- func(x, b, p, a)
    }
    if (fun == "Delta" || fun == "ACO") {
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], a)
        IndTest[, 2] <- func(x, a)
    }
    
    if (fun == "Pxx") {
        if (is.null(d)) 
            d <- distance(spatobj = spatobj, folder = folder, shape = shape)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], d, fdist)
        IndTest[, 2] <- func(x, d, fdist)
    }
    
    if (fun == "ACE") {
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        if (is.null(dc)) 
            dc <- distcenter(spatobj = spatobj, folder = folder, shape = shape, center)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], a, dc)
        IndTest[, 2] <- func(x, a, dc)
    }
    
    if (fun == "HShannon" || fun == "JPielou" || fun == "ISimpson" || 
          fun == "GiniMulti" || fun == "DMulti" || fun == "PIsol" || 
          fun == "RelDivers") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]])
        IndTest[, 2] <- func(x)
    }
    
    
    for (i in 1:nvar) IndTest[i, 3] <- mean(resim[i, ])
    
    for (i in 1:nvar) {
        ncontor <- 0
        for (k in 1:nsim) if (IndTest[i, 2] > resim[i, k]) 
            ncontor <- ncontor + 1
        IndTest[i, 4] <- ncontor + 1
        IndTest[i, 5] <- max((nsim - ncontor)/(nsim + 1), 1/(nsim + 1))
    }
    return(IndTest)
} 
