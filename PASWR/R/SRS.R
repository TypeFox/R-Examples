SRS <-
function(POPvalues, n)
  { # SRS generates all possible SRS's of size n
    # from the population in vector POPvalues
    # by calling the function Combinations.
    N <- length(POPvalues)
    store <- t(Combinations(N, n))
    matrix(POPvalues[t(store)], nrow = nrow(store), byrow = TRUE) }

