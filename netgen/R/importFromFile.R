#' Import a network from proprietary format.
#'
#' @param filename [\code{character(1)}]\cr
#'   File name.
#' @return Nothing
#' @export
importFromFile = function(filename) {
  assertFile(filename, access = "r")

  fh = file(filename, open = "r")
  on.exit(close(fh))

  network = list()
  network = readSpecificationPart(fh, network)

  # convert stuff to numeric
  n.dim = as.numeric(network$dimension)
  n.nodes = as.numeric(network$n_nodes)
  n.depots = as.numeric(network$n_depots)
  n.clusters = as.numeric(network$n_clusters)

  # number of lines to skip
  line.nr = length(network) + 1L

  # read the data block to a data frame
  df = read.table(file = filename, header = TRUE, skip = line.nr, sep = ",")

  coord.ids = paste0("x", 1:n.dim)

  depot.coordinates = NULL
  if (n.depots > 0L) {
    depot.coordinates = as.matrix(df[1:n.depots, coord.ids])
    df = df[-c(1:n.depots), , drop = FALSE]
  }

  network = makeNetwork(
    name = network$name,
    comment = network$comment,
    coordinates = as.matrix(df[, coord.ids]),
    membership = df$membership,
    depot.coordinates = depot.coordinates,
    lower = as.numeric(network$lower),
    upper = as.numeric(network$upper)
  )

  if (!is.null(df$arrival.time)) {
    network$arrival.times = df$arrival.time
  }
  return(network)
}
