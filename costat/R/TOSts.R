TOSts <-
function(spec){
#
# Calculate Test Statistic on spectrum spec
#

#
# Get dimensions
#
J <- spec$nlevels
n <- 2^J

#
# Turn spectrum coefficients into matrix
#
m <- matrix(spec$D, nrow=J, ncol=n, byrow=TRUE)
#
# Compute empirical variance for each level
#
answer <- sum(apply(m,1,var))

return(answer)
}
