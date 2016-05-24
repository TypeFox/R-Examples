#'@title Get neighbours to geohashes
#'@description Geohashes are calculated using a fixed-sized box, which means it's easy
#'to take a geohash and extract the neighbouring boxes around it, in each direction.
#'These functions either extract all neighbours, or individual neighbours, depending on
#'your preference.
#'
#'@aliases neighbours
#'@rdname neighbours
#'
#'@param hashes a vector of geohashes, which can be computed with \code{\link{gh_encode}}
#'
#'@seealso \code{\link{gh_encode}} for generating hashes and \code{\link{gh_decode}} for resolving
#'them into latitude and longitude pairs.
#'
#'@examples
#'#Get a single neighbours
#'north("ezs42")
#'
#'#Get all neighbours!
#'gh_neighbours("ezs42")
#'@export
north <- function(hashes){
  return(gh_neighbour(hashes, c(1, 0)))
}

#'@rdname neighbours
#'@export
northeast <- function(hashes){
  return(gh_neighbour(hashes, c(1, 1)))
}

#'@rdname neighbours
#'@export
east <- function(hashes){
  return(gh_neighbour(hashes, c(0, 1)))
}

#'@rdname neighbours
#'@export
southeast <- function(hashes){
  return(gh_neighbour(hashes, c(-1, 1)))
}

#'@rdname neighbours
#'@export
south <- function(hashes){
  return(gh_neighbour(hashes, c(-1, 0)))
}

#'@rdname neighbours
#'@export
southwest <- function(hashes){
  return(gh_neighbour(hashes, c(-1, -1)))
}

#'@rdname neighbours
#'@export
west <- function(hashes){
  return(gh_neighbour(hashes, c(0, -1)))
}

#'@rdname neighbours
#'@export
northwest <- function(hashes){
  return(gh_neighbour(hashes, c(1, -1)))
}
