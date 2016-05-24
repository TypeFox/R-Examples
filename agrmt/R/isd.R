isd <-
function(V, tolerance=0.1) {
  # Classification over time (ISD-system)
  # Arguments:         V = vector with values at 3 time points (transition points)
  #            tolerance = tolerance (absolute)
  # Example: V <- c(10,10,0)
  if (!length(V) == 3) {warning("Warning: no vector with three time points.") # input validation
      r <- list(type= NA, description = NA) # returning NA for both type and skew (so we can use ISD()$type)
      return(r)}
  x    <- NULL
  x[1] <- compareValues(V[1],V[2], tolerance=tolerance)
  x[2] <- compareValues(V[2],V[3], tolerance=tolerance)
  if (x[1]==1  & x[2]==1)  I <- list(type=1, description="increase, increase")
  if (x[1]==1  & x[2]==0)  I <- list(type=2, description="increase, flat")
  if (x[1]==1  & x[2]==-1) I <- list(type=3, description="increase, decrease")
  if (x[1]==0  & x[2]==1)  I <- list(type=4, description="flat, increase")
  if (x[1]==0  & x[2]==0)  I <- list(type=5, description="flat, flat")
  if (x[1]==0  & x[2]==-1) I <- list(type=6, description="flat, decrease")
  if (x[1]==-1 & x[2]==1)  I <- list(type=7, description="decrease, increase")
  if (x[1]==-1 & x[2]==0)  I <- list(type=8, description="decrease, flat")
  if (x[1]==-1 & x[2]==-1) I <- list(type=9, description="decrease, decrease")
  return(I)
  }
