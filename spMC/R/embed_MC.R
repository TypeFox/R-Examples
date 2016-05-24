embed_MC <-
function(data, coords, loc.id, direction) {
  # Transition probabilities estimation for Embedded Markov Chain 
  #
  #      data vector of data
  #    coords matrix of coordinates
  #    loc.id location Id (which_lines output)
  # direction vector (or versor) of choosen direction

  if (!is.matrix(coords)) coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  nc <- dim(coords)[2]
  if (length(direction) != nc) stop("wrong length of direction vector")
  if (!is.factor(data)) data <- as.factor(data)
  nk <- nlevels(data)
  labels <- levels(data)
  n <- length(data)
  if (n < nk^2) stop("there are not enough data to estimate the parameters")
  if (n != dim(coords)[1]) stop("the number of data is not equal to the number of coordinates")
  if (n != length(loc.id)) stop("length of \"loc.id\" must be equal to the data length")
  loc.id <- as.integer(loc.id)
  
  ord <- order(abs(direction), decreasing = TRUE)
  ord <- cbind(loc.id, coords[, ord])
  ord <- lapply(apply(ord, 2, list), unlist)
  ord <- do.call("order", ord)
  data <- data[ord]
  loc.id <- loc.id[ord]

  tcount <- vector("integer", nk^2)
  tcount <- .C('cEmbedTrans', n = as.integer(n), nk = as.integer(nk), 
               locId = as.integer(loc.id), data = as.integer(data), 
               tcount = as.integer(tcount), PACKAGE = "spMC")$tcount
  storage.mode(tcount) <- "double"
  tcount <- .C('embedTProbs', nk = as.integer(nk), tp = as.double(tcount), 
               PACKAGE = "spMC")$tp
  tcount <- matrix(tcount, ncol = nk)
  diag(tcount) <- NA
  colnames(tcount) <- rownames(tcount) <- labels
  return(tcount)
}

