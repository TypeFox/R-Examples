sorted <- function(x) {

################################################################################
# Function: sorted
# Programmer: Tom Kincaid
# Date: August 13, 2001
# Description:
#   This function determines whether the input set of values is a nondecreasing
#   sequence.
#   Input is a set of values.
#   Output is a logical variable, where TRUE = sorted and FALSE = not sorted.
#   Other Functions Required: None
################################################################################

# Calculate additional required values

   n <- length(x)

# Determine whether the input set of values is a nondecreasing sequence

   sorted <- max(order(x) - (1:n)) == 0

# Return the result

   sorted
}
