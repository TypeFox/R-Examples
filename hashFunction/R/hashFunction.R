#' Calculate hash value of input x using SpookyHash
#' http://en.wikipedia.org/wiki/Jenkins_hash_function
#' @param x a number or a character
#' @return an integer hash value
spooky.32 <- function(x) {
  if (is.integer(x)) {
    storage.mode(x) <- "integer"
    return (.Call("spooky32int", x, PACKAGE ="hashFunction"))
  } else if (is.character(x)){
    return (.Call("spooky32str", x, PACKAGE ="hashFunction"))
  } else {
    warning("Only integer type or character type are supported")
    return (NA);
  }
}

#' Calculate hash value of input x using CityHash
#' https://code.google.com/p/cityhash/
#' @param x a number or a character
#' @return an integer hash value
cityhash.64 <- function(x) {
  if (is.integer(x)) {
    storage.mode(x) <- "integer"
    return (.Call("cityhash64int", x, PACKAGE ="hashFunction"))
  } else if (is.character(x)){
    return (.Call("cityhash64str", x, PACKAGE ="hashFunction"))
  } else {
    warning("Only integer type or character type are supported")
    return (NA);
  }
}

#' Calculate hash value of input x using Murmur3
#' https://code.google.com/p/cityhash/
#' @param x a number or a character
#' @return an integer hash value
murmur3.32 <- function(x) {
  if (is.integer(x)) {
    storage.mode(x) <- "integer"
    return (.Call("murmur3_32int", x, PACKAGE ="hashFunction"))
  } else if (is.character(x)){
    return (.Call("murmur3_32str", x, PACKAGE ="hashFunction"))
  } else {
    warning("Only integer type or character type are supported")
    return (NA);
  }
}


