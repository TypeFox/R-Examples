compareValues <-
function(A,B,tolerance=0.1) {
  # Helper function: How do two values compare?
  # Arguments: A = value 1
  #            B = value 2
  #    tolerance = tolerance (absolute value, default 0.1)
  if (A > B) r <- -1 else r <- 1 # A is bigger: relationship -1; B is bigger: relationship 1
  if (abs(A-B)-tolerance <= 0) r <- 0
  return(r)
  }
