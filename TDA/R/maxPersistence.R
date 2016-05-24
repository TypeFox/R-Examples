maxPersistence <-
function(FUN, parameters, X, lim, by, maxdimension = length(lim) / 2 - 1,
         sublevel = TRUE, library = "Dionysus", B = 30, alpha = 0.05,
         bandFUN = "bootstrapBand", distance = "bottleneck",
         dimension = min(1, maxdimension), p = 1, parallel = FALSE,
         printProgress = FALSE, weight = NULL) {

  if (!is.function(FUN)) {
    stop("FUN should be a function")
  }
  if (!is.numeric(parameters)) {
    stop("parameters should be a numeric vector")
  }
  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(lim) || length(lim) %% 2 != 0) {
    stop("lim should be either a numeric matrix or a numeric vector of even elements")
  }
  if (!is.numeric(by) || any(by <= 0)) {
    stop("by should be positive")
  }
  if (2 * NCOL(X) != length(lim)) {
    stop("dimension of X does not match with lim")
  }
  if (length(by) != 1 && length(by) != NCOL(X)) {
    stop("by should be either a number or a vector of length equals dimension of grid")
  }
  if (!is.numeric(maxdimension) || length(maxdimension) != 1 || maxdimension < 0) {
    stop("maxdimnsion should be a nonnegative integer")
  }
  if (!is.logical(sublevel)) {
    stop("sublevel should be logical")
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library <- "Dionysus"
  }
  if (library == "phat" || library == "Phat") {
    library <- "PHAT"
  }
  if (library != "Dionysus" && library != "PHAT") {
    stop("library should be a string: either 'Dionysus' or 'PHAT'")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B should be a positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    stop("alpha should be a number between 0 and 1")
  }
  if (bandFUN != "bootstrapBand" && bandFUN != "bootstrapDiagram") {
    stop("bandFUN should be a string: either 'bootstrapBand' or 'bootstrapDiagram'")
  }
  if (distance != "wasserstein" && distance != "bottleneck") {
    stop("distance should be a string: either 'bottleneck' or 'wasserstein'")
  }
  if (!is.numeric(dimension) ||
      any(dimension < 0 || dimension > maxdimension)) {
    stop("dimension should be a integer or a vector of integer, with the range between 0 and maxdimension")
  }
  if (!is.numeric(p) || length(p) != 1 || p < 1) {
    stop("p should be a positive integer")
  }
  if (!is.logical(parallel)) {
    stop("parallel should be logical")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }
  if (!is.null(weight) && (!is.numeric(weight) ||
      (length(weight) != 1 && length(weight) != NROW(X)))) {
    stop("weight should be either NULL, a number, or a vector of length equals the number of sample")
  }

  X <- as.matrix(X)
  maxdimension <- min(maxdimension, NCOL(X) - 1)
  Grid <- gridBy(lim = lim, by = by)[["grid"]]
  
  Kseq <- length(parameters)
  eps <- numeric(Kseq)
  numberSignificant <- numeric(Kseq)
  significantPers <- numeric(Kseq)
  Pers <- list()
  
  if (printProgress) {
    cat("0   10   20   30   40   50   60   70   80   90   100\n")
    cat("|----|----|----|----|----|----|----|----|----|----|\n")
    cat("*")
  }
  percentageFloor <- 0

  for (i in seq(along = parameters)){

    if (is.null(weight)) {
      Diag <- gridDiag(X = X, FUN = FUN, lim = lim, by = by, FUNvalues = NULL,
          maxdimension = maxdimension, sublevel = sublevel, library = library,
          location = FALSE, printProgress = FALSE, diagLimit = NULL,
          parameters[i])[["diagram"]]
    } else {
      Diag <- gridDiag(X = X, FUN = FUN, lim = lim, by = by, FUNvalues = NULL,
          maxdimension = maxdimension, sublevel = sublevel, library = library,
          location = FALSE, printProgress = FALSE, diagLimit = NULL,
          weight = weight, parameters[i])[["diagram"]]
    }
    Diag[1,3] <- Diag[1,2] #remove first component with infinite persistence
    Pers[[i]] <- cbind(Diag[,1], Diag[,3] - Diag[,2])
    colnames(Pers[[i]]) <- c("dimension", "Persistence")
      
    if (bandFUN == "bootstrapDiagram") {
      eps[i] <- bootstrapDiagram(X = X, FUN = FUN, lim = lim, by = by,
          maxdimension = maxdimension, sublevel = sublevel, library = library,
          B = B, alpha = alpha, distance = distance, dimension = dimension,
          p = p, parallel = parallel, printProgress = FALSE, weight = weight,
          parameters[i])
      selctDim <- rep(FALSE, NROW(Pers[[i]]))
      for (dim in dimension) {
        selctDim <- selctDim | (Pers[[i]][,1] == dim)
      }
      numberSignificant[i] <- sum(Pers[[i]][selctDim, 2] > (2 * eps[i]))
      significantPers[i] <- sum(pmax(0, Pers[[i]][selctDim, 2] - (2 * eps[i])))

    } else {
      eps[i] <- bootstrapBand(X = X, FUN = FUN, Grid = Grid, B = B,
          alpha = alpha, parallel = parallel, printProgress = FALSE,
          weight = weight, parameters[i])[["width"]]
      numberSignificant[i] <- sum(Pers[[i]][, 2] > (2 * eps[i]))
      significantPers[i] <- sum(pmax(0, Pers[[i]][, 2] - (2 * eps[i])))
    }
    
    if (printProgress) {
      for (j in seq_len(floor((100 * i / Kseq - percentageFloor) / 2))) {
        cat("*")
        percentageFloor <- percentageFloor + 2
      }
    }
  }
  if (printProgress)
  { 
    cat("\n")
  }

  # Two criterions
  Param1 <- parameters[which(numberSignificant == max(numberSignificant))]
  Param2 <- parameters[which(significantPers == max(significantPers))]
  
  out <- list("parameters" = parameters, "sigNumber" = numberSignificant,
      "sigPersistence" = significantPers, "bands" = eps, "Persistence" = Pers,
      "bandFUN" = bandFUN, "dimension" = dimension)
  class(out) <- "maxPersistence"
  attributes(out)[["call"]] <- match.call()
  return (out)
}
