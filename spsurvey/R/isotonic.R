isotonic <- function(y, minval, maxval) {

################################################################################
# Function: isotonic
# Programmer: Tom Kincaid
# Date: August 15, 2001
# Description:
#   This function performs isotonic regression of an input set of values so that
#   the output set of values is a nondecreasing sequence.  The output set of
#   values is truncated to the range: [minval, maxval]. 
#   Input:
#      y = the set of values on which to perform isotonic regression.
#      minval = the minimum value for the output set of values.
#      maxval = the maximum value for the output set of values.
#   Output is the revised set of values
#   Other Functions Required:
#      sorted - determine whether a set of values is a nondecreasing
#         sequence
################################################################################

# Calculate additional required values

   n <- length(y)

# Truncate the output set of values to the range: [minval, maxval]

   y <- pmax(y, minval)
   y <- pmin(y, maxval)

# Perform isotonic regression

   while (!sorted(y)) {
      i <- 1
      while (i < n) {
         j <- i
         while (y[j] >= y[j+1] & j < n)
            j <- j+1
         if (i < j) y[i:j] <- mean(y[i:j])
         i <- j+1
      }
   }

# Return the result

   y
}
