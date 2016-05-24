getlen <-
function(data, coords, loc.id, direction, zero.allowed = FALSE) {
  # Empirical estimation of lengths (for embeded data)
  #
  #         data vector of data
  #       coords matrix of coordinates
  #       loc.id location Id (which_lines output)
  #    direction vector (or versor) of choosen direction
  # zero.allowed logical, it allows zero stratum lengths

  if (!is.matrix(coords)) coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  nc <- dim(coords)[2]
  if (length(direction) != nc) stop("wrong length of direction vector")
  if (!is.factor(data)) data <- as.factor(data)
  nk <- nlevels(data)
  n <- length(data)
  if (n != dim(coords)[1]) stop("the number of data is not equal to the number of coordinates")
  if (n != length(loc.id)) stop("length of \"loc.id\" must be equal to the data length")
  storage.mode(loc.id) <- "integer"

  ord <- order(abs(direction), decreasing = TRUE)
  ord <- cbind(loc.id, coords[, ord])
  ord <- lapply(apply(ord, 2, list), unlist)
  ord <- do.call("order", ord)
  data <- data[ord]
  loc.id <- loc.id[ord]
  coords <- coords[ord, ]
  lens <- .C('cEmbedLen', n = as.integer(n), nc = as.integer(nc),
             coords = as.double(coords), locId = as.integer(loc.id),
             data = as.integer(data), cemoc = as.integer(vector("integer", n)),
             maxcens = as.double(vector("numeric", n)), tlen = as.double(vector("numeric", n)),
             PACKAGE = "spMC")[c(1, 6:8)]
  mycenslen <- lens$maxcens[1:lens$n]
  lens$maxcens <- NULL
  lens$categories <- as.factor(levels(data)[lens$cemoc[1:lens$n]])
  lens$length <- lens$tlen[1:lens$n]
  lens$maxcens <- mycenslen
  if (!zero.allowed) {
    idx <- lens$length != 0
    lens$categories <- lens$categories[idx]
    lens$length <- lens$length[idx]
    lens$maxcens <- lens$maxcens[idx]
  }
#   else {
#     lens$length <- lens$length + lens$maxcens
#   }
  lens$direction <- direction
  lens$zeros <- zero.allowed
  lens$cemoc <- NULL
  lens$tlen <- NULL
  lens$n <- NULL
  class(lens) <- "lengths"
  return(lens)
}
