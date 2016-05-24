patternAgreement <-
function(P, old=FALSE) {
   # Calculates agreement A from a pattern vector
   # Arguments:   P = pattern vector
   #            old = use old Unimodality algorithm (passed on)
   # Example: Second pattern on p.333: P <- c(1,0,1,0,0,0,0)
   if (max(P) > 1) stop("Error: Input is not a pattern vector (only 0 and 1 are allowed).")
   # This error should never occur unless the function is called directly.
   K <- length(P)                     # number of categories
   # Counting triplets
   TDU <- 0                           # begin with count = 0, prepare
   TU  <- 0                           # begin with count = 0, prepare
   for (i in 1:(K-2)) {               # repeat for position A
      for (j in (i+1):(K-1)) {        # repeat for position B
         for (m in (j+1):K) {         # repeat for position C
	    if (P[i] == 1 & P[j] == 0 & P[m] == 1) TDU <- TDU + 1 # 101 pattern, bimodal (TDU)
	    if (P[i] == 1 & P[j] == 1 & P[m] == 0) TU <- TU + 1   # 110 pattern, unimodal (TU)
	    if (P[i] == 0 & P[j] == 1 & P[m] == 1) TU <- TU + 1   # 011 pattern, unimodal (TU)
	    # all other patterns are not counted
	 }
      }
   }
   if (old==TRUE) (U <- (TU-TDU)/(TU+TDU))           # using the old algorithm (outlined in endnotes)
   else (U <- ((K-2)*TU-(K-1)*TDU)/((K-2)*(TU+TDU))) # normal case: U as in equation (2) on p.332
   S <- sum(P)                        # number of non-empty
   A <- U*(1-(S-1)/(K-1))             # calculating agreement A
   if (is.nan(A)) (A <- 0)            # lack of agreement, defined as 0
   if (sum(P) == 1) (A <- 1)          # only one value, defined as 1
   return(A)
   }
