##' Compute all pairwise distances for a population. This function
##' is simply a wrapper for \code{dist} that returns only a matrix
##'
##' @title Distance Matrix for a mating scene
##' @param scene a matingScene object
##' @return a matrix of all pairwise comparisons with attributes for order of
##' identifiers (idOrder)
##' @seealso \code{\link{dist}}
##' @examples
##' pop <- simulateScene()
##' distance <- pairDist(pop)
pairDist <- function(scene) {
  distMat <- as.matrix(dist(scene[, c("x", "y")]))
  attr(distMat, "idOrder") <- scene$id
  attr(distMat, "dimnames") <- NULL
  distMat
}

##' Find the k nearest neighbors for all individuals in a population. This
##' function is simply a wrapper for \code{FNN::knn.dist}.
##'
##' @title Get k Nearest Neighbors
##' @param scene a matingScene object
##' @param k integer of how many nearest neighbors to get
##' @return a matrix where the rows are all individuals and the columns are
##' their k nearest neighbors
##' @export
##' @seealso \code{\link{knn.dist}}
##' @examples
##' pop <- simulateScene(10)
##' kNearNeighbors(pop, 3)
kNearNeighbors <- function(scene, k) {
  knnMatrix <- FNN::knn.dist(scene[c("x", "y")], k = k, algorithm = "brute")
  rownames(knnMatrix) <- scene$id
  colnames(knnMatrix) <- paste("k", 1:k, sep = "")
  knnMatrix
}

##' Calculate one of several measures of spatial proximity
##'
##' @title Make potentials object--spatial proximity
##' @param scene a matingScene object
##' @param method one of "maxProp", and "maxPropSqrd" see details for
##' further description
##' @param proximityFun a function used to calculate proximity. Not yet
##' implemented
##' @param averageType whether to calculate individual and population proximity
##' using the mean or median
##' @param subject whether you want pair, individual, population, or all.
##' Specifying more than one is allowed.
##' @param zeroPotDist the distance at which potential should be equal to zero
##' @return A potentials object containing one more more of the following, depending the
##' input for \code{subject}: \cr
##' If \code{subject} is "population" the return list will contain a numeric
##' value that has a range depending on the \code{method}. If
##' \code{subject} is "pair" the return list will contain a matrix
##' with all pairwise proximity comparisons. If \code{subject} is "individual"
##' the return list will contain a dataframe with a column containing IDs and
##' a column containing proximity values. If \code{subject} is "all"
##' the return list will contain all three of the items above.
##' @details If \code{method} is "maxProp" then proximity between two
##' individuals will be calculated as 1 - distance/max(distance).
##' If \code{method} is "maxPropSqrd" then proximity between two
##' individuals will be calculated as (1 - distance/max(distance))^2.
##' @author Danny Hanson
##' @examples
##' pop <- simulateScene()
##' proximity(pop, "maxProp")
proximity <- function(scene, method, proximityFun = NULL, averageType = "mean",
                      subject = "all", zeroPotDist = NULL) {
  method <- match.arg(method, c("maxProp", "maxPropSqrd"))
  subject <- match.arg(subject, c("all", "pair", "population", "individual"),
                       several.ok = TRUE)
  if (is.list(scene) & !is.data.frame(scene)) {
    potential <- lapply(scene, proximity, method, proximityFun, averageType, subject)
  } else {
    n <- nrow(scene)
    distMatrix <- pairDist(scene)
    if (averageType == "mean") {
      average <- mean
    } else if (averageType == "median") {
      average <- median
    }
    if (is.null(zeroPotDist)) {
      zeroPotDist <- max(distMatrix)
    }
    # deal with pop size
    if (n < 2) {
      stop("Can't calculate proximity for population size less than 2")
    }

    if (method == "maxProp") {
      if (subject %in% c("all", "pair")) {
        pairProx <- 1 - distMatrix/zeroPotDist
      }

      distNoDiag <- distMatrix[-seq(1, n^2, n+1)]
      distNoDiag[distNoDiag > zeroPotDist] <- zeroPotDist
      distNoDiagMat <- matrix(distNoDiag, nrow = n, byrow = T)
      pairProx2 <- 1 - distNoDiagMat/zeroPotDist

      indProx <- data.frame(id = scene$id, proximity = -1)
      if (averageType == "mean") {
        indProx$proximity <- rowMeans(pairProx2)
      } else if (averageType == "median") {
        indProx$proximity <- row_medians(pairProx2)
      }

      popProx <- average(indProx[,2])
    } else if (method == "maxPropSqrd") {
      if (subject %in% c("all", "pair")) {
        pairProx <- (1 - distMatrix/zeroPotDist)^2
      }

      distNoDiag <- distMatrix[-seq(1, n^2, n+1)]
      distNoDiag[distNoDiag > zeroPotDist] <- zeroPotDist
      distNoDiagMat <- matrix(distNoDiag, nrow = n, byrow = T)
      pairProx2 <- (1 - distNoDiagMat/zeroPotDist)^2

      indProx <- data.frame(id = scene$id, proximity = -1)
      if (averageType == "mean") {
        indProx$proximity <- rowMeans(pairProx2)
      } else if (averageType == "median") {
        indProx$proximity <- row_medians(pairProx2)
      }

      popProx <- average(indProx[,2])
    }

    # return
    potential <- list()
    if ("population" %in% subject) {
      potential$pop <- popProx
    }
    if ("individual" %in% subject) {
      potential$ind <- indProx
    }
    if ("pair" %in% subject) {
      potential$pair <- pairProx
    }
    if ("all" %in% subject) {
      potential$pop <- popProx
      potential$ind <- indProx
      potential$pair <- pairProx
    }
    attr(potential, "t") <- FALSE
    attr(potential, "s") <- TRUE
    attr(potential, "c") <- FALSE
    potential
  }
}
