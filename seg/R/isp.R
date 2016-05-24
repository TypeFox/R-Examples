# ------------------------------------------------------------------------------
# isp()
# ------------------------------------------------------------------------------
isp <- function(x, data, nb, fun, verbose = FALSE, ...) {

  # The input object 'x' can be either a class of 'Spatial' (or its associates)
  # or 'data.frame'. Depending on the class of 'x', 'data' may not be required.
  # The internal function 'chksegdata()' processes the information as required
  # by the current function.
  if (verbose)
    tmp <- chksegdata(x, data)
  else
    tmp <- suppressMessages(chksegdata(x, data))
  coords <- tmp$coords; pdf <- tmp$data;
  rm(tmp)

  # Verify 'coords' and 'data'
  if (ncol(pdf) < 2)
    stop("'data' must be a matrix with at least two columns", call. = FALSE)
  else if (!is.numeric(pdf))
    stop("'data' must be a numeric matrix", call. = FALSE)
  else if (nrow(pdf) != nrow(coords))
    stop("'data' must have the same number of rows as 'x'", call. = FALSE)

  # If the user did not specify 'dist', do the followings:
  if (missing(nb))
    nb <- as.matrix(dist(x = coords, ...))
  
  # If 'fun' is not given, use the default (i.e., negative exponential):
  if (missing(fun))
    fun <- function(z) { exp(-z) }
  
  # If the distance matrix is symmetric and if all the left-to-right diagonals
  # are 0, we can save some time:
  if (isSymmetric(nb)) {
    pairID <- t(combn(1:nrow(nb), 2))             # Pairs of census tracts
    pairDist <- as.numeric(as.dist(nb))           # Distance between pairs
  } else { # If not symmetric:
    pairID <- expand.grid(1:nrow(nb), 1:nrow(nb)) # Pairs of census tracts
    pairDist <- as.numeric(nb)                    # Distance between pairs
  }
  VALID <- which(pairDist != 0)
  pairID <- pairID[VALID,]
  pairDist <- pairDist[VALID]
  
  # Determines the degree of spatial interaction based on distances
  speffect <- fun(pairDist)
  
  pRow <- apply(pdf, 1, sum) # Total population by census tracts
  pCol <- apply(pdf, 2, sum) # Total population by groups
  
  pA <- list()
  pB <- sum(pRow[pairID[,1]] * pRow[pairID[,2]] * speffect) / sum(pCol)
  for (i in 1:ncol(pdf))
    pA[[i]] <- sum(pdf[pairID[,1],i] * pdf[pairID[,2],i] * speffect) / pCol[i]
  
  as.vector(sum(unlist(pA)) / pB)
}
