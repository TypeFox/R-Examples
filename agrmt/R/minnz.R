minnz <-
function(V) {
   # Calculates the smallest value of the vector except for 0 (non-zero minumum)
   # Argument: vector
   C <- NULL        # prepare
   k <- length(V)   # count to
   for (i in 1:k) { # ceck all
      if ((V[i] == 0) == FALSE) (C[i] <- V[i]) else (C[i] <- 9999919) # if V[i] is not 0, add it to C
   }
   m <- min(C)               # minimum of V, not counting 0
   if (max(V) == 1) (m <- 1) # fix for binary vectors (0,1)
   if (m == 9999919) (warning("Error: Minimum calculation failed."))  # warning because of hard-coded replacement
   return(m)
   }
