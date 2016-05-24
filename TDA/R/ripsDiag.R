ripsDiag <-
function(X, maxdimension, maxscale, dist = "euclidean", library = "GUDHI",
         location = FALSE, printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(maxdimension) ||
      length(maxdimension) != 1 || maxdimension < 0) {
    stop("maxdimnsion should be a nonnegative integer")
  }
  if (!is.numeric(maxscale) || length(maxscale) != 1) {
    stop("maxscale should be a number")
  }
  if (dist != "euclidean" && dist != "arbitrary") {
    stop ("dist should be either 'euclidean' or 'arbitrary'")
  }
  if (library == "gudhi" || library == "Gudhi") {
    library <- "GUDHI"
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library <- "Dionysus"
  }
  if (library == "phat" || library == "Phat") {
    library <- "PHAT"
  }
  if (library != "GUDHI" && library != "Dionysus" && library != "PHAT") {
    stop("library should be a string: either 'GUDHI, 'Dionysus', or 'PHAT'")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }

  if (dist == "arbitrary" && library == "GUDHI") {
    library <- "Dionysus"
  }

  X <- as.matrix(X)

  max_num_pairs <- 5000  # to be added as an option

  # in 32bit architectures Dionysus/PHAT L2 doesn't work
  #ripsOut <- RipsDiag(X = X, maxdimension = maxdimension,
  #    maxscale = maxscale, dist = dist, library = library,
  #    location = location, printProgress = printProgress)
  if (dist == "euclidean" && library != "GUDHI" &&
      .Machine[["sizeof.pointer"]] != 8) {
    ripsOut <- RipsDiag(X = as.matrix(dist(X)), maxdimension = maxdimension,
        maxscale = maxscale, dist = "arbitrary", library = library,
        location = location, printProgress = printProgress)
  } else {
    ripsOut <- RipsDiag(X = X, maxdimension = maxdimension,
        maxscale = maxscale, dist = dist, library = library,
        location = location, printProgress = printProgress)
  }

  if (location == TRUE) {
    if (dist == "euclidean") {
      BirthLocation <- X[ripsOut[[2]][, 1], ]
      DeathLocation <- X[ripsOut[[2]][, 2], ]
      if (library == "Dionysus")
      {
        CycleLocation <- lapply(ripsOut[[3]], function(c) {X[c, ]})
      }
    } else {
      BirthLocation <- ripsOut[[2]][, 1]
      DeathLocation <- ripsOut[[2]][, 2]
      if (library == "Dionysus")
      {
        CycleLocation <- ripsOut[[3]]
      }
    }
  }

  Diag <- ripsOut[[1]]
  if (NROW(Diag) > 0) {
    ## change Inf values to maxscale
    Diag[which(Diag[, 3] == Inf), 3] <- maxscale  
  }

  colnames(Diag) <- c("dimension", "Birth", "Death")
  class(Diag) <- "diagram"
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  attributes(Diag)[["scale"]] <- c(0, maxscale)
  attributes(Diag)[["call"]] <- match.call()
  if (location == FALSE || library == "GUDHI")
  {
    out <- list("diagram" = Diag)
  } else if (library == "PHAT")
  {
    out <- list("diagram" = Diag, "birthLocation" = BirthLocation,
        "deathLocation" = DeathLocation)
  } else
  {
    out <- list("diagram" = Diag, "birthLocation" = BirthLocation,
        "deathLocation" = DeathLocation, "cycleLocation" = CycleLocation)
  }
  return (out)
}
