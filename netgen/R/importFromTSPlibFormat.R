#' Import network from (extended) TSPlib format.
#'
#' @note The extended TSPlib contains additional specification parts and a cluster
#' membership section. Currently only the import of symmetric TSP instances is possible.
#'
#' @param filename [\code{character(1)}]\cr
#'   Path to TSPlib file.
#' @param round.distances [\code{logical(1)}]\cr
#'   Should the distances of EUC_2D instances be rounded to the nearest integer value?
#'   Default is \code{TRUE}.
#' @return [\code{Network}]
#'   Network object.
#' @export
importFromTSPlibFormat = function(filename, round.distances = TRUE) {
  requirePackages("stringr", why = "netgen::importFromTSPlibFormat")
  assertFile(filename, access = "r")
  assertFlag(round.distances)

  fh = file(filename, open = "r")
  on.exit(close(fh))

  network = list()
  network = readSpecificationPart(fh, network)
  network = deepCheckSpecification(network)

  n.points = as.integer(network$dimension)

  line = str_trim(readLines(fh, 1L))
  while (length(line) > 0 && line != "EOF" && line != "" && !is.na(line)) {
    if (line == "NODE_COORD_SECTION") {
      network[["coordinates"]] = readNodeCoordinates(fh, n.points)
    }
    if (line == "DISPLAY_DATA_SECTION") {
      network[["display_data"]] = readNodeCoordinates(fh, n.points)
    }
    if (line == "CLUSTER_MEMBERSHIP_SECTION") {
      network[["membership"]] = readClusterSection(fh, n.points)
    }
    if (line == "EDGE_WEIGHT_SECTION") {
      network = readEdgeWeightsSection(fh, network, n.points)
    }
    line = str_trim(readLines(fh, 1L))
  }

  # postprocessing
  network$edge_weights = getNetworkEdgeWeights(network, round.distances)
  network$coordinates = getNetworkCoordinates(network)

  # finally generate netgen {Clustered}Network object
  makeNetwork(
    name = network$name,
    comment = network$comment,
    coordinates = network$coordinates,
    distance.matrix = network$edge_weights,
    lower = if (!is.null(network$lower)) as.numeric(network$lower) else NULL,
    upper = if (!is.null(network$upper)) as.numeric(network$upper) else NULL,
    membership = network$membership,
    edge.weight.type = network$edge_weight_type
  )
}

deepCheckSpecification = function(network) {
  assertList(network)
  network$dimension = as.integer(network$dimension)
  if (is.null(network$dimension)) {
    stopf("TSPlib format error: Mandatory DIMENSION specification is missing.")
  }
  if (!is.null(network$lower)) {
    network$lower = as.numeric(network$lower)
  }
  if (!is.null(network$upper)) {
    network$upper = as.numeric(network$upper)
  }
  if (!is.character(network$name) || network$name == "") {
    stopf("NAME is mandatory and connot be empty.")
  }
  if (network$type != "TSP") {
    stopf("At the moment only the symmetric TSP files are supported, but you provided TYPE '%s'", network$type)
  }
  return(network)
}

# Extract specifications.
#
# @note Multiple comments possible.
#
# @param fh [connection]
#   File handle.
# @param network [list]
#   List of key-values pairs.
# @return [list]
#   Modified network.
readSpecificationPart = function(fh, network) {
  repeat {
    line = readLines(fh, 1L)
    line.parts = strsplit(line, split = "[[:space:]]*:[[:space:]]*")[[1]]

    # we reached the SECTIONs part
    if (length(line.parts) != 2L) {
      pushBack(line, fh)
      break
    }
    key = tolower(line.parts[1])
    value = str_trim(line.parts[2])
    # multiple comments are allowed. We store all of them and not just
    # the last one
    if (key == "comment" && !is.null(network[[key]])) {
      network[[key]] = c(network[[key]], value)
    } else {
      network[[key]] = value
    }
  }
  return(network)
}

# Construct final coordinates.
#
# @param network [list]
#   List of key-values pairs.
# @return [list]
#   Modified network.
getNetworkCoordinates = function(network) {
  # if NODE_COORD_SECTION was available
  if (!is.null(network$coordinates)) {
    return(network$coordinates)
  }
  # if NO_DISPLAY is set
  if (!is.null(network$display_data_type)) {
    if (network$display_data_type == "NO_DISPLAY") {
      stopf("There are no coordinates available for the instance '%s'.", network$name)
    }
  }
  # if DISPLAY_DATA section was provided
  if (!is.null(network$display_data)) {
    return(network$display_data)
  }

  # otherwise try to reconstruct cooridinates based on distance matrix
  edge_weights = network$edge_weights
  if (!is.null(edge_weights)) {
    return(cmdscale(dist(edge_weights), k = 2))
  }
  stopf("No coordinates available and no possibility to guess them for the given instance '%s'.", network$name)
}

# Construct final distance matrix.
#
# @param network [list]
#   List of key-values pairs.
# @return [list]
#   Modified network.
getNetworkEdgeWeights = function(network, round.distances) {
  edge_weights = network$edge_weights
  # if there is an EDGE_WEIGHT_SECTION
  if (!is.null(edge_weights)) {
    return(edge_weights)
  }
  ewt = network$edge_weight_type
  coordinates = network$coordinates
  n.points = as.integer(network$dimension)
  mapping = list("EUC_2D" = "euclidean", "CEIL_2D" = "euclidean", "MAX_2D" = "maximum", "MAN_2D" = "manhattan")
  # if there is a nice EDGE_WEIGHT_TYPE
  if (ewt %in% names(mapping)) {
    distance.matrix = as.matrix(dist(coordinates, method = mapping[[ewt]]))
    if (ewt != "CEIL_2D" && round.distances) {
      distance.matrix = round(distance.matrix)
    }
    # occasionally we need to ceil
    if (ewt == "CEIL_2D") {
      # round up to the next integer
      distance.matrix = ceiling(distance.matrix)
    }
    return(distance.matrix)

  } else if (ewt == "ATT") {
    # special "pseudo-Euclidean" distance (as it is called in Reinelt 95)
    distance.matrix = matrix(0, ncol = n.points, nrow = n.points)
    for (i in 1:n.points) {
      for (j in 1:n.points) {
        if (i == j) {
          next
        }
        # See Reinelt 95
        distance.matrix[i, j] = sqrt(sum((coordinates[i, ] - coordinates[j, ])^2) / 10.0)
        tmp = getNextInteger(distance.matrix[i, j])
        if (tmp < distance.matrix[i, j]) {
          distance.matrix[i, j] = tmp + 1L
        } else {
          distance.matrix[i, j] = tmp
        }
      }
    }
    return(distance.matrix)
  } else if (ewt == "GEO") {
    # geographical distance
    coordinates = network$coordinates
    x = coordinates[, 1]
    degrees = floor(x)
    min = x - degrees
    latitude = pi * (degrees + 5 * min / 3) / 180
    y = coordinates[, 2]
    degrees = floor(y)
    min = y - degrees
    longitude = pi * (degrees + 5 * min / 3) / 180
    earth.radius = 6378.3888
    distance.matrix = matrix(0, ncol = n.points, nrow = n.points)
    for (i in 1:n.points) {
      for (j in 1:n.points) {
        if (i == j) {
          next
        }
        q1 = cos(longitude[i] - longitude[j])
        q2 = cos(latitude[i] - latitude[j])
        q3 = cos(latitude[i] + latitude[j])
        distance.matrix[i, j] = round(earth.radius * acos(0.5 * ((1 + q1) * q2 - (1 - q1) * q3)) + 1)
      }
    }
    return(distance.matrix)
  }
  stopf("Unsupported EDGE_WEIGHT_TYPE: '%s'.", ewt)
}

# "Round" an numeric value to its next integer value.
# E.g. getNextInteger(1.4) == 1 but getNextInteger(1.6) == 2.
#
# @param x [numeric]
#   Numeric vector.
# @return [integer]
#   Integer vector.
getNextInteger = function(x) {
  as.integer(floor(x + 0.5))
}

# Helper function.
#
# Takes a file handle and the number of nodes and
# returns the "raw" matrix of coordinates.
#
# @param fh [connection]
#   File handle.
# @param n [integer(1)]
#   Problem dimension, i.e., the number of nodes.
# @return [matrix]
#   (n x 2) matrix of coordinates.
readNodeCoordinates = function(fh, n) {
  # <integer> <real> <real>
  raw.coordinates = scan(fh, nmax = 3 * n, quiet = TRUE)
  # get rid of node number (every third element)
  raw.coordinates = raw.coordinates[-seq(1, 3 * n, by = 3)]
  coordinates = matrix(raw.coordinates, ncol = 2L, byrow = TRUE)
  return(coordinates)
}

# Helper function.
#
# Takes a file handle, a network and the number of nodes and
# returns the "raw" distance matrix.
#
# @param fh [connection]
#   File handle.
# @param network [list]
#   List of key-value pairs.
# @param n [integer(1)]
#   Problem dimension, i.e., the number of nodes.
# @return [matrix]
#   (n x n) distance matrix.
readEdgeWeightsSection = function(fh, network, n) {
  ewt = network$edge_weight_type
  ewf = network$edge_weight_format
  if (is.null(ewt)) {
    stopf("Edge weight section found, but not edge weight type specified.")
  }
  if (ewt != "EXPLICIT") {
    stopf("Currently only explicit edge weight types are supported.")
  }
  if (is.null(ewf)) {
    stopf("Edge weight section is found, but no edge weight format given.")
  } else if (ewf == "FULL_MATRIX") {
    #FIXME: all the read function have the same signature
    # Construct a mapping from EDGE_WEIGHT_TYPE to the corresponding
    # function to make the code nicer
    edge.weights = readExplicitEdgeWeights(fh, n)
  } else if (ewf == "UPPER_ROW") {
    edge.weights = readUpperRowWeights(fh, n)
  } else if (ewf == "LOWER_ROW") {
    edge.weights = readLowerRowWeights(fh, n)
  } else if (ewf == "UPPER_DIAG_ROW") {
    edge.weights = readUpperDiagRowWeights(fh, n)
  } else if (ewf == "LOWER_DIAG_ROW") {
    edge.weights = readLowerDiagRowWeights(fh, n)
  } else {
    #FIXME: add support for the remaining types
    # UPPER_COL, LOWER_COL, UPPER_DIAG_COL and LOWER_DIAG_COL
    stopf("Unsupported EDGE_WEIGHT_FORMAT: '%s'.", ewf)
  }
  network$edge_weights = edge.weights
  return(network)
}

# Reads distance matrix for EDGE_WEIGHT_FORMAT == UPPER_DIAG_ROW
#
# Takes a file handle, the number of nodes and
# construct the "raw" distance matrix from an upper row represenation
# with diagonal elements.
#
# @param fh [connection]
#   File handle.
# @param n [integer(1)]
#   Problem dimension, i.e., the number of nodes.
# @return [matrix]
#   (n x n) distance matrix.
readUpperDiagRowWeights = function(fh, n) {
  distance.matrix = matrix(0, ncol = n, nrow = n)
  distances = scan(fh, nmax = (n * (n + 1)) / 2, quiet = TRUE)
  i = 1L
  j = 1L
  for (k in 1:length(distances)) {
    distance.matrix[i, j] = distances[k]
    j = j + 1L
    if (j > n) {
      i = i + 1L
      j = i
    }
  }
  distance.matrix[lower.tri(distance.matrix)] = t(distance.matrix)[lower.tri(distance.matrix)]
  return(distance.matrix)
}

# Reads distance matrix for EDGE_WEIGHT_FORMAT == LOWER_DIAG_ROW
#
# Takes a file handle, the number of nodes and
# construct the "raw" distance matrix from an lower row represenation
# with diagonal elements.
#
# @param fh [connection]
#   File handle.
# @param n [integer(1)]
#   Problem dimension, i.e., the number of nodes.
# @return [matrix]
#   (n x n) distance matrix.
readLowerDiagRowWeights = function(fh, n) {
  distance.matrix = matrix(0, ncol = n, nrow = n)
  distances = scan(fh, nmax = (n * (n + 1)) / 2, quiet = TRUE)
  i = 1L
  j = 1L
  for (k in 1:length(distances)) {
    distance.matrix[i, j] = distances[k]
    j = j + 1L
    if (j > i) {
      i = i + 1L
      j = 1L
    }
  }
  distance.matrix[upper.tri(distance.matrix)] = t(distance.matrix)[upper.tri(distance.matrix)]
  return(distance.matrix)
}

# Reads distance matrix for EDGE_WEIGHT_FORMAT == UPPER_ROW
#
# Takes a file handle, the number of nodes and
# construct the "raw" distance matrix from an upper row represenation.
#
# @param fh [connection]
#   File handle.
# @param n [integer(1)]
#   Problem dimension, i.e., the number of nodes.
# @return [matrix]
#   (n x n) distance matrix.
readUpperRowWeights = function(fh, n) {
  distance.matrix = matrix(0, ncol = n, nrow = n)
  distances = scan(fh, nmax = (n * (n - 1)) / 2, quiet = TRUE)
  i = 1L
  j = 2L
  for (k in 1:length(distances)) {
    distance.matrix[i, j] = distances[k]
    j = j + 1L
    if (j > n) {
      i = i + 1L
      j = i + 1L
    }
  }
  distance.matrix[lower.tri(distance.matrix)] = t(distance.matrix)[lower.tri(distance.matrix)]
  return(distance.matrix)
}

# Reads distance matrix for EDGE_WEIGHT_FORMAT == LOWER_ROW
#
# Takes a file handle, the number of nodes and
# construct the "raw" distance matrix from an lower row represenation.
#
# @param fh [connection]
#   File handle.
# @param n [integer(1)]
#   Problem dimension, i.e., the number of nodes.
# @return [matrix]
#   (n x n) distance matrix.
readLowerRowWeights = function(fh, n) {
  distance.matrix = matrix(0, ncol = n, nrow = n)
  distances = scan(fh, nmax = (n * (n - 1)) / 2, quiet = TRUE)
  i = 2L
  j = 1L
  for (k in 1:length(distances)) {
    distance.matrix[i, j] = distances[k]
    j = j + 1L
    if (j > i) {
      i = i + 1L
      j = 1L
    }
  }
  distance.matrix[upper.tri(distance.matrix)] = t(distance.matrix)[upper.tri(distance.matrix)]
  return(distance.matrix)
}

# Reads distance matrix for EDGE_WEIGHT_FORMAT == EXPLICIT
#
# Takes a file handle, the number of nodes and
# construct the "raw" distance matrix from an upper diagonal represenation.
#
# @param fh [connection]
#   File handle.
# @param n [integer(1)]
#   Problem dimension, i.e., the number of nodes.
# @return [matrix]
#   (n x n) distance matrix.
readExplicitEdgeWeights = function(fh, n) {
  distances = scan(fh, nmax = n * n, quiet = TRUE)
  matrix(distances, ncol = n, nrow = n, byrow = TRUE)
}

# Reads membership vector.
#
# Takes a file handle, the number of nodes and
# construct the vector of cluster membership.
#
# @param fh [connection]
#   File handle.
# @param n [integer(1)]
#   Problem dimension, i.e., the number of nodes.
# @return [integer]
#   Integer vector of cluster membership.
readClusterSection = function(fh, n) {
  membership = as.integer(scan(fh, nmax = n, quiet = TRUE))
  return(membership)
}
