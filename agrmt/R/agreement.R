agreement <-
function(V, old=FALSE) {
   # Calculates agreement A for multiple layers
   # Arguments:  V = frequency vector
   #           old = use old Unimodality algorithm
   # Example: V <- c(30,40,210,130,530,50,10)
   if (length(V) < 3 ) {
      warning("Warning: length of vector < 3, agreement A is not defined.")
      return(NA)}
   if (min(V) < 0) stop("Error: negative values found in frequency vector.")
   # This error only occurs if the input is not a frequency vector. Use collapse() to generate a frequency vector.
   AA <- 0        # begin with empty agreement A (overall), prepare
   k <- length(V) # number of categories
   N <- sum(V)    # number of cases
   R <- V         # remainder
   for (i in 1:k) {                  # repeat for each layer i
      P <- patternVector(R)          # get the pattern vector for layer i
      if (max(P) == 0) break         # remainder is empty, all layers are analyzed
      A <- patternAgreement(P, old)  # agreement A for layer i
      m <- minnz(R)                  # get non-zero minimum of remainder R
      L <- P * m                     # layer i with the values
      w <- sum(L)/N                  # weight of layer i
      AA <- AA + w * A               # adding agreement of layer i to overall agreement
      R <- R - L                     # new reminder
   }
   return(AA)
   }
