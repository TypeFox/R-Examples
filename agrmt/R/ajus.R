ajus <-
function(V, tolerance=0.1, variant="modified") {
  # Determine the shape of a distribution
  # Arguments:         V = frequency vector
  #            tolerance = tolerance (absolute value)
  #              variant = "strict" for AJUS, "modified" for AJUSFL
  # Example: V <- c(30,40,210,130,530,50,10)
  if (min(V) < 0) stop("Error: negative values found in frequency vector.") # input validation
  # (1) identify patterns: 0 flat, 1 increase, -1 decrease
  n <- length(V)    # number of items
  z <- n-1          # number of breaks
  if (n < 3) {warning("Warning: too few values to classify distribution.") # input validation
      r <- list(type=NA, skew=NA) # returning NA for both type and skew (so we can use AJUS()$type)
      return(r)}
  x <- NULL         # prepare
  for (i in 1:z) {
    x[i] <- compareValues(V[i], V[i+1], tolerance=tolerance)
    }
  # Example vector V gives x = c(1  1 -1  1 -1 -1)
  # or with tolerance=.5:  x = c(0  1  0  1 -1 -1)
  # (2) identify shape
  min.x <- min(x)
  max.x <- max(x)
  if (min.x == 0 & max.x == 0 & compareValues(V[1], V[n], tolerance=tolerance) == 0) A <- "F" else {    # flat distribtion, type F not in AJUS
    if (max.x <  1) A <- "J" else {      # no 1,  thus only 0 or + 1 (single peak at left end)  -- "L"
      if (min.x > -1) A <- "J"  else {   # no -1, thus only 0 or 1   (single peak at right end) -- "J"
        xs <- reduceVector(x)            # remove 0 and repeated values; not use unique(), because I want same values at different positions ("type S")
        if (isTRUE(all.equal(xs, c(1,-1)))) A <- "A" else {
          # isTRUE(all.equal(V[i],m)
          if (isTRUE(all.equal(xs, c(-1,1)))) A <- "U" else A <- "S"
        }
      }
    }
  }
  # F if: flat, no peak;                polarization; all 0
  # A if: unimodal, peak in the middle; consensus;    0 or +1, then 0 or -1
  # J if: unimodal, peak at either end; consensus;    0 or + 1, xor 0 or -1
  # U if: bimodal,  peak at both ends;  polarization; 0 or -1, then 0 or +1
  # S if: bimodal,  multiple peaks;     polarization; else
  # (3) identify skew
  m <- round(n/2,0) # midpoint of vector V
  S <- compareValues(sum(V[1:m]),sum(V[m:n]), tolerance=tolerance) # S = skew [-1,0,+1] (negative, symmetric, positive)
  # manually set skew for "J" type distributions (necessary because long tails can change value S):
	if (A == "J" & max.x < 1)  S <- -1 # single peak at left end  -- "L"
	if (A == "J" & min.x > -1) S <- 1  # single peak at right end -- "J"
  # "strict" versus "modified"
  if (variant == "strict" & A == "F") A <- NA              # "F" not defined in AJUS strict
  if (variant == "modified" & A == "J" & S == -1) A <- "L" # "L" as modified type
  r <- list(type = A, skew = S, pattern = x)
  return(r)
  }
