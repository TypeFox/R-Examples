# created: Lutz Prechelt,  1997-10-13
# changed: Lutz Prechelt,  2006-03-01

#-------------------- data handling functions:

eqc <- function(x, number=4, overlap=1/sqrt(length(x))) {
  equal.count(x,number,overlap)
}

grep.b <- function(pattern, x) {
  # like grep(), with two differences:
  #  - returns a boolean match vector, not a vector of match indices.
  #  - coerces x to character (important for factors)
  # This is useful for 'subset' expressions in Lattice calls.
  result <- rep(FALSE, times=length(x))
  result[grep(pattern, as.character(x))] <- TRUE
  return (result)
}

grep.s <- function(pattern, x) {
  grep(pattern, x, value=TRUE)
}

orderavg <- function(x, by) {
  # orders the observations from vector 'x' using the ranking criteria given
  # in the matrix 'by' (length(x) rows, n colums), a 'by' vector is coerced.
  # The ranks computed from each matrix column are summed row-wise,
  # the vector is sorted according to these rank sums, and returned.
  # This function can be used to pair elements from two vectors (by ranking
  # them both by equivalent criteria).
  # TODO: add weighting of criteria
  rnk.x    <- apply(cbind(by),    2, rank) # rank each criterion
  rnksum.x <- apply(rbind(rnk.x), 1, sum)  # combine criteria
  return(x[order(rnksum.x)])        # re-order criteria
}

tableNA <- function(...) table(..., exclude = NULL)

or.else <- function(x, alternative = 0) {
  # returns x where available and alternative where is.na(x)
  ifelse(is.na(x), alternative, x)
}

#-------------------- debug helper

printn <- function(x, digits) {
  # for debugging: prints out the name of the first argument and its value.
  # will bomb if the digits argument is used and x cannot be
  # coerced into numeric
  if(!missing(digits))
    x <- round(x, digits)
  print(match.call()[[2]])
  print(x)
  invisible()
}

tracebck = function() {
  # Only one line per call, as this suppresses large arguments.
  tb = get(".Traceback", envir = baseenv())
  sapply(tb, function(x) x[[1]])
}
