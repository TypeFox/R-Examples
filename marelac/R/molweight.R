## -----------------------------------------------------------------------------
## Molar Weight; converts from mol to g
## -----------------------------------------------------------------------------

molweight <- function(species) {
  if (!is.vector(species)) stop("species must be a vector")
  if (!is.character(unlist(species))) stop("species must be character")
  molweight <- function(species) {
    # atomicweight is the 'named list-version' of the IUPAC data
    with(marelac::atomicweight, {
      ## insert * before number (with one or more digits)
      s1 <- gsub("([0-9]+)", "*\\1+", species)
      ## insert + after capital letters
      s1 <- gsub("([a-z,A-Z])", "\\1+", s1, perl = TRUE)
      ## remove + before lower case letters
      s1 <- gsub("\\+([a-z])", "\\1", s1, perl = TRUE)
      ## replace +* with only *
      s1 <- gsub("+*", "*", s1, fixed = TRUE)
      ## remove trailing +
      s1 <- gsub("\\+$", "", s1)
      ## remove trailing +)
      s1 <- gsub("\\+\\)", ")", s1)
      ## calculate molar mass
      eval(parse(text = s1))
    })
  }
  sapply(species, molweight)
}

