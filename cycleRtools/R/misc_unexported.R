# These functions aren't exported, but are otherwise useful.
#  ------------------------------------------------------------------------

# Useful binary operator for subsetting.
"%btwn%" <- function(x, rng) x >= rng[[1]] & x <= rng[[2]]

# As opposed to %in%.
"%notin%" <- function(x, table) match(x, table, nomatch = 0L) == 0L

# For partial name matching.
"%grep_in%" <- function(x, table) any(grepl(x, table))

# Generate colour palette without dependencies.
base_pal <- function(x, s) {
  # s = unique values;
  # x = subsequently matched against s to give hex colour codes.
  if (length(s) == 1)
    "#000000"
  else
    rainbow(n = length(s))[match(x, s)]
}

# Error message for functions with specific cycleRdata methods.
format_error <- function() {
  message(paste(
    "Data is not properly formatted.",
    "Import the data via a read_ride function with the 'format' argument as TRUE.",
    "Returning NULL.", sep = "\n"))
  return(NULL)
}

# Interpolate NAs; NB: won't replace NAs at the start of input.
na_approx <- function(x) {
  approx(x    = seq_len(length(x))[!is.na(x)],
         y    = x[!is.na(x)],
         xout = seq_len(length(x)))$y
}

# Unload DynLib.
.onUnload <- function (libpath) {
  library.dynam.unload("cycleRtools", libpath)
}
