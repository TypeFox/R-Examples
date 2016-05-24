#' Read in a TSPLIB style Traveling Salesman Problem from a
#' file.
#'
#' The current state of the parser does not understand all variants
#' of the TSPLIB format. Much effort has been spent making the
#' parser as robust as possible. It will stop as soon as it sees
#' input it cannot handle.
#'
#' @param path [\code{character(1)}]\cr
#'   Character string containing path to file in TSPLIB format.
#' @return [\code{\link{tsp_instance}}].
#' @export
read_tsplib_instance = function(path) {
  assertFile(path, access = "r")

  con = file(path, open = "r")
  on.exit(close(con))

  tsp_instance = list()

  # Read specification part
  tsp_instance = read_tsplib_specification(con, tsp_instance)

  if (substr(tsp_instance$TYPE, 1, 3) != "TSP")
    stop("Currently the only supported TYPE is: TSP!")

  # Read data part
  while (next_line_exists(con)) {
    line = next_line(con)

    ## What section is this?
    if (line == "NODE_COORD_SECTION") {
      tsp_instance = read_tsplib_node_coords(con, tsp_instance)
    } else if (line == "EDGE_WEIGHT_SECTION") {
      tsp_instance = read_tsplib_edge_weights(con, tsp_instance)
    } else if (line == "DISPLAY_DATA_SECTION") {
      tsp_instance = read_tsplib_display_data(con, tsp_instance)
    } else if (line == "FIXED_EDGES_SECTION") {
      tsp_instance = read_tsplib_fixed_edges(con, tsp_instance)
    } else if (line == "") {
      next
    } else if (line == "EOF") {
      break
    } else {
      stop("Unhandled data section '", line, "' in TSPLIB-file.")
    }
  }
  # wrap in object of type \code{tsp_instance}.
  tsp_instance$EDGE_WEIGHTS = edge_weights(tsp_instance)
  coords = node_coordinates(tsp_instance)
  tsp_instance(coords = as.matrix(coords), dists = tsp_instance$EDGE_WEIGHTS)
}

# Read next line from a connection.
#
# @param con [\code{connection}]\cr
#   Connection object.
# @return [\code{character(1)}] Next line of text from \code{con} or \code{NULL} if no
#   more lines of text are available.
next_line = function(con) {
  str_trim(readLines(con, n = 1, warn = FALSE))
}

# Peek at next line of text from a connection.
#
# @param con [\code{connection}]\cr
#   Connection object.
# @return [\code{character(1)}] Next line of text from \code{con} or \code{NULL} if no
#   more lines of text are available. The line is not removed
#   from the connection so the next call to \code{next_line} will
#   return it.
peek_line = function(con) {
  line = readLines(con, 1, warn = FALSE)
  pushBack(line, con)
  line
}

# Check if one more line exists.
#
# @param con [\code{connection}]\cr
#   Connection object.
# @return [\code{logical(1)}] Boolean value indicating
#   whether another line in connection exists.
next_line_exists = function(con) {
  length(peek_line(con)) > 0
}


read_tsplib_specification = function(con, tsp_instance) {
  line = next_line(con)
  # All entries in specification part are of the form
  # <keyword> : <value>
  while (str_detect(line, ":")) {
    specification = unlist(strsplit(line, "[[:space:]]*:[[:space:]]*"))
    tsp_instance[[specification[1]]] = specification[2]
    line = next_line(con)
  }
  tsp_instance$DIMENSION = as.integer(tsp_instance$DIMENSION)

  # Last line read does _not_ belong to tsp_instance, so push it
  # back onto the connection.
  pushBack(line, con)
  tsp_instance
}

read_tsplib_node_coords = function(con, tsp_instance) {
  if (is.null(tsp_instance$NODE_COORD_TYPE)) {
    ## W00t - breakin the law!
    ##
    ## Seriously: We blindly guess the type of coords that follow.
    tsp_instance$NODE_COORD_TYPE = "TWOD_COORDS"
  }

  if (tsp_instance$NODE_COORD_TYPE == "TWOD_COORDS") {
    if (is.null(tsp_instance$DIMENSION))
      stop("Unknown dimension of TSP instance in function 'read_tsplib_node_coords'.")
    # read all the coordinates
    node_coords = read.table(con, header = FALSE, nrows = tsp_instance$DIMENSION)
    ## remove first col, which just contains index numbers 1:n
    tsp_instance$NODE_COORDS = node_coords[,-1]
    tsp_instance$DIMENSION = nrow(node_coords)
  } else {
    stop("Unsupported NODE_COORD_TYPE '", tsp_instance$NODE_COORD_TYPE, "'.")
  }
  tsp_instance
}

read_tsplib_edge_weights = function(con, tsp_instance) {
  if (is.null(tsp_instance$EDGE_WEIGHT_TYPE))
    stop("EDGE_WEIGHT_TYPE not specified but EDGE_WEIGHT_SECTION present!")

  dimension = tsp_instance$DIMENSION
  if (tsp_instance$EDGE_WEIGHT_TYPE == "EXPLICIT") {
    if (tsp_instance$EDGE_WEIGHT_FORMAT == "LOWER_DIAG_ROW") {
      number_of_edge_weights = dimension * (dimension + 1) / 2
      edge_weights = scan(con, n = number_of_edge_weights, quiet = TRUE)
      # Delete diagonal
      edge_weights = edge_weights[-cumsum(1:dimension)]

      # Build symmetric distance matrix
      m = matrix(0, nrow = dimension, ncol = dimension)
      m[upper.tri(m)] = edge_weights
      edge_weights = t(m)[lower.tri(m)]

      attr(edge_weights, "Diag") = FALSE
      attr(edge_weights, "Upper") = FALSE
    } else if (tsp_instance$EDGE_WEIGHT_FORMAT == "UPPER_DIAG_ROW") {
      number_of_edge_weights = dimension * (dimension + 1) / 2
      edge_weights = scan(con, n = number_of_edge_weights, quiet = TRUE)
      edge_weights = edge_weights[-cumsum(c(1, rev(2:dimension)))]
      attr(edge_weights, "Diag") = FALSE
      attr(edge_weights, "Upper") = TRUE
    } else if (tsp_instance$EDGE_WEIGHT_FORMAT == "UPPER_ROW") {
      number_of_edge_weights = dimension * (dimension - 1) / 2
      edge_weights = scan(con, n=number_of_edge_weights, quiet = TRUE)
      attr(edge_weights, "Diag") = FALSE
      attr(edge_weights, "Upper") = TRUE
    } else if (tsp_instance$EDGE_WEIGHT_FORMAT == "FULL_MATRIX") {
      number_of_edge_weights = tsp_instance$DIMENSION ** 2
      edge_weights = scan(con, n = number_of_edge_weights, quiet = TRUE)
      tsp_instance$edge_weights = edge_weights
    } else if (tsp_instance$EDGE_WEIGHT_FORMAT == "FUNCTION") {
      stop("Encountered 'function' edge weight format which is currently not supported.")
    } else {
      stop("Unsupported EDGE_WEIGHT_FORMAT '", tsp_instance$EDGE_WEIGHT_FORMAT, "'.")
    }
    if (tsp_instance$EDGE_WEIGHT_FORMAT == "FULL_MATRIX") {
      edge_weights = as.dist(matrix(edge_weights, dimension))
    } else {
      attr(edge_weights, "Size") = dimension
      class(edge_weights) = "dist"
    }
    tsp_instance$EDGE_WEIGHTS = edge_weights
  } else {
    stop("Unsupported EDGE_WEIGHT_TYPE '", tsp_instance$EDGE_WEIGHT_TYPE, "'.")
  }
  tsp_instance
}

read_tsplib_display_data = function(con, tsp_instance) {
  if (is.null(tsp_instance$DISPLAY_DATA_TYPE )) {
    ## W00t - breakin the law!
    ##
    ## Seriously: Fixup missing DISPLAY_DATA_TYPE field:
    tsp_instance$DISPLAY_DATA_TYPE <- "TWOD_DISPLAY"
  }

  if (tsp_instance$DISPLAY_DATA_TYPE == "TWOD_DISPLAY") {
    dimension = tsp_instance$DIMENSION
    display_data = read.table(con, header=FALSE, nrows=tsp_instance$DIMENSION)
    ## remove first col, which just contains index numbers 1:n
    tsp_instance$DISPLAY_DATA = display_data[, -1]
  }
  tsp_instance
}

read_tsplib_fixed_edges = function(con, tsp_instance) {
  ## Ignore for now:
  i = scan(con, n = 1, quiet = TRUE)
  while (i > 0) {
    i = scan(con, n = 1, quiet = TRUE)
  }
  tsp_instance
}

# Return the edge weights (distances) of a TSP
#
# @param tsp_instance instance of a TSP.
#
# @return [\code{dist}] Object containing the distances.
#
# @references G. Reinelt TSPLIB 95.
#
edge_weights = function(tsp_instance) {
  edge_weight_type = tsp_instance$EDGE_WEIGHT_TYPE
  if (!is.null(tsp_instance$EDGE_WEIGHTS)) {
    tsp_instance$EDGE_WEIGHTS
  } else if (edge_weight_type == "EUC_2D") {
    dist(node_coordinates(tsp_instance))
  } else if (edge_weight_type == "CEIL_2D") {
    ceiling(dist(node_coordinates(tsp_instance)))
  } else if (edge_weight_type == "ATT") {
    ## We cheat here and ignore the scaling factor:
    ceiling(dist(node_coordinates(tsp_instance)))
  } else if (edge_weight_type == "GEO") {
    # See tsplib95 format explanation for details on this.
    degrees_to_geo = function(degrees) {
      full_degrees = floor(degrees)
      minutes = degrees - full_degrees
      pi * (full_degrees + 5 * minutes / 30) / 180
    }
    coordinates = node_coordinates(tsp_instance)
    latitude = degrees_to_geo(coordinates[,1])
    longitude = degrees_to_geo(coordinates[,1])
    ## FIXME: Inefficient!
    distance_matrix = outer(1:length(latitude), 1:length(latitude),
                             function(i, j) {
                               earth_radius = 6378.3888
                               q1 = cos(longitude[i] - longitude[j])
                               q2 = cos(latitude[i] - latitude[j])
                               q3 = cos(latitude[i] + latitude[j])
                               round(earth_radius * acos(0.5 * ((1 + q1) * q2 - (1 - q1) * q3)) + 1)
                             })
    edge_weights = distance_matrix[upper.tri(distance_matrix)]
    attr(edge_weights, "Diag") = FALSE
    attr(edge_weights, "Upper") = TRUE
    attr(edge_weights, "Size") = length(latitude)
    class(edge_weights) = "dist"
    edge_weights
  } else {
    stop("Unsupported EDGE_WEIGHT_TYPE '", edge_weight_type, "'.")
  }
}

# Return (2D) coordinates of each node
#
# @param tsp_instance instance of a TSP.
#
# @return A matrix with two columns (x and y coordinate) and as many
# rows as there are nodes in the TSP.
node_coordinates = function(tsp_instance) {
  if (!is.null(tsp_instance$NODE_COORDS)) {
    tsp_instance$NODE_COORDS
  } else if (!is.null(tsp_instance$DISPLAY_DATA_TYPE)
             && tsp_instance$DISPLAY_DATA_TYPE == "NO_DISPLAY") {
    stop("No node coordinates available (DISPLAY_DATA_TYPE = 'NO_DISPLAY').")
  } else if (!is.null(tsp_instance$DISPLAY_DATA)) {
    tsp_instance$DISPLAY_DATA
  } else if (!is.null(tsp_instance$EDGE_WEIGHTS)) {
    edge_weights = tsp_instance$EDGE_WEIGHTS
    if (any(edge_weights == 0)) {
      ## Guess coordinates since we only have a distance matrix
      warning("Some nodes have distance 0.")
      cmdscale(edge_weights, k = 2)
    } else {
      projection = sammon(edge_weights, niter = 1000, tol = 1e-9, trace = FALSE)
      ##FIXME: The subsequent warning is converted into an error (gr17.tsp)
      if (projection$stress > 0.0001) {
        ## We cannot find a 2D embedding which reporduces the distance
        ## matrix to decent accuracy. Too bad.
        warning("Node coordinates are only approximate (stress=",
                projection$stress, ")")
      }
      projection$points
    }
  } else {
    stop("No node coordinates available.")
  }
}

#' Read in a TSPLIB style Traveling Salesman Problem tour from a file
#'
#' @param path [\code{character(1)}]\cr
#'   Filename of file containing a TSP tour.
#' @return [\code{\link[TSP]{TOUR}}]
#'   TOUR object from package TSP, containing order of cities, tour length and
#'   method name that generated this solution.
read_tsplib_tour = function(path) {
  # determine file extension
  splitted = strsplit(path, ".", fixed = TRUE)[[1]]
  file_ext = splitted[length(splitted)]
  con = file(path, open = "r")
  on.exit(close(con))

  tsp_tour = list()
  # tour is given in tsplib format
  if(file_ext == "tour") {
    dimension = NULL

    line = next_line(con)
    while (str_trim(peek_line(con)) != "TOUR_SECTION") {
        field = str_split(line, ":")[[1]]
        field_name = str_trim(field[1])
        field_value = str_trim(field[2])
        # save dimension (neccessary to extract tour )
        if(field_name == "DIMENSION") {
          dimension = as.integer(field_value)
        }
        line = next_line(con)
    }

    while (next_line_exists(con)) {
      line = next_line(con)

      # what section is it?
      if (line == "TOUR_SECTION") {
        tsp_tour = read_tsplib_tour_section(con)
      }
    }

  # tour is given in concorde format
  } else if(file_ext == "sol") {
    line = next_line(con)
    dimension = as.integer(line)
    tsp_tour = numeric()
    while (next_line_exists(con)) {
      line = next_line(con)
      sub_tour = as.integer(unlist(strsplit(line, " ")[[1]]))
      tsp_tour = c(tsp_tour, sub_tour)
    }
    # concorde returns nodes 0 to (n-1), but we need nodes between 1 and n
    tsp_tour = tsp_tour + 1
  } else {
    stop("BAM! Unknown tour file format.")
  }
  return(as.integer(tsp_tour))
}

# read_tsplib_tour_section - extracts the nodes (cities) on the TSP tour
#
# @param con - open connection to file
# @param dimension
#
# @return vector of nodes on the TSP tour.
read_tsplib_tour_section = function(con) {
  nodes = numeric()
  line = next_line(con)
  while(line != "-1") {
    # we have to do some creepy stuff here to achieve correct tour
    temp = as.integer(c(str_split(str_trim(str_split(line, "  ")[[1]]), " "),
                        recursive = TRUE))
    nodes = c(nodes, temp)
    line = next_line(con)
  }
  return(nodes)
}
