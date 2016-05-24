################################################################### 
# THIS FILE CONTAINS FUNCTIONS THAT THAT IMPLEMENT THE "LOF"      #
# ALGORITHM (Breunig et. al. 2000)                                #
#                                                                 #
# ** IMPORTANT NOTICE **                                          #
# THIS CODE IS ALMOST A 100% COPY OF THE CODE PREVIOUSLY AVAILABLE#
# IN PACKAGE dprep (Acuna et. al. 2009) THAT WAS REDRAWN FROM THE #
# CRAN REPOSITORY.                                                #
###################################################################
# Author : Luis Torgo (ltorgo@inescporto.pt) BUT SEE NOTICE ABOVE #
# License: GPL                                                    #
###################################################################

# =================================================================
# REFERENCES:
# Acuna, E., and Members of the CASTLE group at UPR-Mayaguez, (2009).
#   dprep: Data preprocessing and visualization functions for classification. R
#   package version 2.1.
#
# Breunig, M., Kriegel, H., Ng, R., and Sander, J. (2000). LOF: identifying
#   density-based local outliers. In ACM Int. Conf. on Management of Data,
#   pages 93–104.
# =================================================================



# =====================================================
# Function that computes the LOF factor for a data set
# using k neighbours
# =====================================================
# Example run:
# data(iris)
# l <- lofactor(iris[,-5],10)
#
lofactor <- function(data,k) {

  data <- as.matrix(data)

  # obtain the k nearest neighbors and their distance from each observation
  distdata <- dist.to.knn(data,k)
  p <- dim(distdata)[2L]
 
  # calculate the local reachability density for each observation in data
  lrddata <- reachability(distdata,k)

  lof <- rep(0,p)

  # computer the local outlier factor of each observation in data
  for (i in 1:p) {
    nneigh <- distdata[2,i]-distdata[1,i]+1
    j <- seq(0,(nneigh-1))
    local.factor <- sum(lrddata[distdata[3+j,i]]/lrddata[i])/nneigh
    lof[i] <- local.factor
  }

  # return lof, a vector with the local outlier factor of each observation
  lof
}




# =====================================================
# Function that returns an object in which each column
# contains the indices of the first k neighours followed
# by the distances to each of these neighbours
# =====================================================
dist.to.knn <- function(dataset,neighbors) {

  numrow <- dim(dataset)[1L]
  mxNN <- neighbors*2+2
  knndist <- matrix(0,nrow=mxNN,ncol=numrow)


  for (i in 1:numrow) {
    # find obervations that make up the k-distance neighborhood for observation dataset[i,]
    neighdist <- knneigh.vect(dataset[i,],dataset,neighbors)
    
    x <- length(neighdist)
    if (x > mxNN) {
      knndist <- rbind(knndist,matrix(rep(0,(x-mxNN)*numrow),ncol=numrow))
      mxNN <- x
    }
    knndist[1:x,i] <- neighdist
  }

  return(knndist[1:mxNN,])
}



# =====================================================
# Function that returns the distance from a vector "x" to   
# its k-nearest-neighbors in the matrix "data"
# =====================================================
knneigh.vect <- function(x,data,k) {

  temp <- as.matrix(data)
  numrow <- dim(data)[1L]
  dimnames(temp) <- NULL
  
  # subtract rowvector x from each row of data
  difference<- scale(temp, x, FALSE)

  # square and add all differences and then take the square root
  dtemp <- drop(difference^2 %*% rep(1, ncol(data)))
  dtemp <- sqrt(dtemp)

  # order the distances
  order.dist <- order(dtemp)
  nndist <- dtemp[order.dist]

  # find distance to k-nearest neighbor
  # uses k+1 since first distance in vector is a 0
  knndist <- nndist[k+1]

  # find neighborhood
  # eliminate first row of zeros from neighborhood 
  neighborhood <- drop(nndist[nndist<=knndist])
  neighborhood <- neighborhood[-1]
  numneigh <- length(neighborhood)

  # find indexes of each neighbor in the neighborhood
  index.neigh <- order.dist[1:numneigh+1]

  # this will become the index of the distance to first neighbor
  num1 <- numneigh+3

  # this will become the index of the distance to last neighbor
  num2 <- numneigh+numneigh+2


  return(c(num1,num2,index.neigh,neighborhood))
}


# =====================================================
# Function that calculates the local reachability density
# of Breuing(2000) for each observation in a matrix, using
# a matrix (distdata) of k nearest neighbors computed by the function dist.to.knn2
# =====================================================
# Example run:
# data(iris)
# 
#
reachability <- function(distdata,k) {

  p <- dim(distdata)[2]
  lrd <- rep(0,p)

  for (i in 1:p) {
    j <- seq(3,3+(distdata[2,i]-distdata[1,i]))
    # compare the k-distance from each observation to its kth neighbor
    # to the actual distance between each observation and its neighbors
    numneigh <- distdata[2,i]-distdata[1,i]+1
    temp <- rbind(diag(distdata[distdata[2,distdata[j,i]],distdata[j,i]]),distdata[j+numneigh,i])

    #calculate reachability
    reach <- 1/(sum(apply(temp,2,max))/numneigh)
    lrd[i] <- reach
  }
  lrd
}

