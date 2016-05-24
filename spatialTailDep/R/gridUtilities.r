# Function checking for integer
is.whole <- function(x, tol = .Machine$double.eps^0.5){abs(x - round(x)) < tol}

regularGrid <- function(height, width) {
  if(! is.whole(height) || ! is.whole(width) || length(height) > 1 || length(width) > 1){
    warning("Height and width must be integers")
  } else{
    locations <- matrix(data = 0, nrow = height * width, ncol = 2)
    locations[, 1] <- rep(1:height, width)
    locations[, 2] <- rep(1:width, each = height)
    return( locations )
  }
}

#' Function \code{selectPairIndices}
#' 
#' Out of a list of locations given by their Cartesian coordinates, selects pairs of locations with a distance as a criterion.
#' 
#' @param locations A \eqn{d} x 2 matrix containing the Cartesian coordinates of \eqn{d} points in the plane.
#' @param maxDistance Pairs of locations with distance not larger than \code{maxDistance} will be selected. If \code{NULL}, then the selection is made based on \code{numberOfPairs}.
#' @param numberOfPairs The number of pairs that will be selected, taking distance as a criterion, with closer pairs being selected. It will return more pairs if there are several equidistant locations.
#' @return A matrix with \code{q} rows and 2 columns, where \code{q} denotes the number of pairs selected. Each row contains the indices of the corresponding pair of selected locations in the input argument \code{locations}.
#' @seealso \code{\link{pairCoordinates}}
#' @export
#' @examples
#' (locations<-cbind(rep(1:4,5),rep(1:5,each=4)))
#' selectPairIndices(locations, maxDistance = 1.5)
selectPairIndices <- function(locations, maxDistance = NULL, numberOfPairs = NULL) {
  if(is.null(maxDistance) && is.null(numberOfPairs)){
    warning("Supply either maxDistance or numberOfPairs")
  } else{
    nloc <- dim(locations)[1]
    ncandidatePairs <- nloc * (nloc - 1) / 2
    candidatePairs <- matrix(0, nrow = ncandidatePairs, ncol = 3)
    index <- 0
    for (i in 1:(nloc-1)) {
      for (j in (i+1):nloc) {
        index <- index + 1
        candidatePairs[index, 1] <- i
        candidatePairs[index, 2] <- j
        candidatePairs[index, 3] <- sqrt(sum((locations[i, ] - locations[j, ])^2))
      }
    }
    if (is.null(maxDistance)) {
      maxDistance <- sort(candidatePairs[, 3])[numberOfPairs]  
    } 
    indices <- candidatePairs[, 3] <= maxDistance
    pairIndices <- candidatePairs[indices, 1:2]
    return(pairIndices)
  }
}

#' Function \code{pairCoordinates}
#' 
#' Given a matrix of coordinates of locations and a matrix of indices of pairs of locations, returns a 
#' matrix with the coordinates of the pairs of locations.
#' 
#' @param locations A \eqn{d} x 2 matrix containing the Cartesian coordinates of \eqn{d} points in the plane.
#' @param pairIndices A \eqn{q} x 2 matrix containing the indices of \eqn{q} pairs of points.
#' @return A \eqn{q} x 4 matrix, where each row gives the Cartesian coordinates of the two locations in the corresponding pair.
#' @seealso \code{\link{selectPairIndices}}
#' @export
#' @examples
#' (locations<-cbind(rep(1:2,3),rep(1:3,each=2)))
#' (pairs <- selectPairIndices(locations, maxDistance = 1.5))
#' pairCoordinates(locations, pairs)
pairCoordinates <- function(locations, pairIndices) {
  q <- nrow(pairIndices)
  result <- matrix(0, nrow = q, ncol = 4)
  for (m in 1:q) {
    result[m, 1:2] <- locations[pairIndices[m, 1], ]
    result[m, 3:4] <- locations[pairIndices[m, 2], ]
  }
  return(result)
}

#' Wind speeds in the Netherlands.
#'
#' Daily maximal speeds of wind gusts, measured in 0.1 m/s. The data are observed at 
#' 22 inland weather stations in the Netherlands. Only the summer months are presented
#' here (June, July, August). Also included are the Euclidian coordinates of the 22 
#' weather stations, where a distance of 1 corresponds to 100 kilometers. For more 
#' information on this dataset, see Einmahl et al. (2014).
#'
#' @format KNMIdata$data is a matrix with 672 rows and 22 columns, KNMI$loc is a matrix with 22 rows
#' and 2 columns.
#' @source \url{http://www.knmi.nl/climatology/daily_data/selection.cgi}
#' @name KNMIdata
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A. and Segers, J. (2014), "An M-estimator of spatial tail dependence". See \url{http://arxiv.org/abs/1403.1975}.
#' @examples
#' data(KNMIdata)
#' locations <- KNMIdata$loc
#' pairIndices <- selectPairIndices(locations, maxDistance = 0.5)
#' Mestimator(KNMIdata$data, locations, pairIndices, k = 60, model="BR",
#' iterate=FALSE, covMat = FALSE)$theta
NULL