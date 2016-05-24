# Force coordinates out of bounds to bounds.
#
# @param coordinates [matrix | data.frame]
#   Coordinates.
# @param out.of.bounds.handling, [character(1)]
#   Strategy to handle out of bounds coordinates.
# @param lower [numeric(1)]
#   Lower box constaint for cube.
# @param upper [numeric(1)]
#   Upper box constaint for cube.
# @return [data.frame]
forceToBounds = function(coordinates, out.of.bounds.handling, lower = 0, upper = 1) {
  if (out.of.bounds.handling == "reset") {
    return(pmin(pmax(coordinates, lower), upper))
  } else if (out.of.bounds.handling == "mirror") {
    for (i in seq(ncol(coordinates))) {
      idx.lower = which(coordinates[, i] < lower)
      coordinates[idx.lower, i] = lower + abs(coordinates[idx.lower, i] - lower)
      idx.upper = which(coordinates[, i] > upper)
      coordinates[idx.upper, i] = upper - abs(coordinates[idx.upper, i] - upper)
    }
    return(coordinates)
  }
  stop("Unknown out.of.bounds.handling!")
}
