patternVector <-
function(V) {
   # Creates pattern vector from frequency vector
   # Argument: V = frequency vector
   k <- length (V) # number of categories
   P <- NULL       # prepare empty P
   for (i in 1:k) {
      if (V[i] == 0) (P[i] <- 0)     # 0 remains 0
      else (P[i] <- 1)               # everything else becomes 1
   }
   return (P)
   }
