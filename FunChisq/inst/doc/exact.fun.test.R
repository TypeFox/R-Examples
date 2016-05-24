## ------------------------------------------------------------------------
require(FunChisq)
# x is a contingency table with row variable for p53 mutation and
#   column variable for CIMP
x <- matrix(c(12,26,18,0,8,12), nrow=2, ncol=3, byrow=TRUE)
x
# Test the functional dependency: p53 mutation -> CIMP
fun.chisq.test(x, method="exact")

