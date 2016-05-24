modes <-
function(V, pos=FALSE, tolerance=0.1) {
  # Find (multiple) modes of frequency vector (roughly the same frequencies)
  # Arguments:    V = frequency vector
  #             pos = positions of categories
  #       tolerance = tolerance to consider values the same (absolute values)
  # Example 1: modes(c(30,40,210,130,530,50,10)) # will find position 5
  # Example 2: modes(c(3,0,4,1))                 # will find position 3
  # Example 3: modes(c(3,0,4,1),pos=-1:2)        # will still find position 3, but give the value of 1 as mode
  # Example 4: modes(c(30,40,500,130,530,50,10),tolerance=30) # will find positions 3 and 5 (500 and 530 nearly same as in the difference smaller or equal to the tolerance)
  if (min(V) < 0) stop("Error: negative values found in frequency vector.")      # input validation
  # This error only occurs if the input is not a frequency vector. Use collapse() to generate a frequency vector.
  if (!pos[1]==FALSE) if (is.numeric(pos) == TRUE) { # input validation: if pos argument is provided
    if (!length(pos) == length(V)) { # pos vector of different length
       warning("Note: length of position vector different from length of frequency vector. Position vector ignored")
       pos <- FALSE } # simply ignore it, better than stopping outright
       }
  k <- length(V)  # number of values
  p <- NULL       # empty for preparation
  if (pos[1]==FALSE) if (is.numeric(pos) == FALSE) pos <- 1:k # no positions provided, assume 1:k (no zero)
  m <- which.max(V) # position of mode
  for (i in 1:k) {  # check each position to see if it is equal to the mode
    # because of the tolerance, frequencies need not be exactly the same
    if (compareValues(V[i],V[m],tolerance=tolerance)==0) p <- c(p,i) # add position of the mode
  }
  # Check whether modes are contiguous (= agreement) (c = {"contiguous", "divided"}
  if ((max(p) - min(p)) < length(p)) c <- TRUE else c <- FALSE
  r <- list(at = pos, frequencies = V, mode = pos[p], positions = p, contiguous = c)
  # at = positions or categories of vector, either given (argument "pos")
  #      or set to 1:k
  # frequencies = frequency vector V (provided)
  # mode = mode(s) of V, considering the categories of the vector
  # positions = position of V at which mode is found
  # contiguous = TRUE/FALSE if modes are contiguous
  return(r)
  }
