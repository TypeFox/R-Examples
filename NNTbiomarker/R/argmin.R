#'  argmin Argmin function for a vector.
#'
#'  Return the index minimizing distance from v to target.
#'  @param v The vector to compare to target.
#'  @param target The value sought in the vector; default=0.
#'  @return  The index in v of the value which is closest to target.

argmin = function(v, target=0)
  sapply(target, function(target)which(abs(v-target) == min(abs(v-target))[1]))
