isMissingValue = function(x) {
  return(is.vector(x) && length(x) == 1L && is.na(x))
}

isMissingName = function(x) {
  return(identical(x, NA_character_))
}
