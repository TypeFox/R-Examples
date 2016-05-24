qj2i <- function(q,j,qmax){
#
# The index "i" is in position (q,j) of an array where the
# first column of the array is 1 to qmax, the second column is
# NA, (qmax + 1) to (2*qmax - 1), the third column is
# NA, NA, (2*qmax):(3*qmax - 3), and so on.
#
if(any(j>q))
    stop("Index \"j\" must be less than or equal to \"q\".\n")
(j-1)*(qmax - j/2) + q
}
