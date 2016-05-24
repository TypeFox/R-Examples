aux3 <- function(m,ind)
{
#
# Auxilliary function to re-structure the information matrix
# in the case where the variances of all components are constrained
# to be equal.
#
x <- apply(m[,ind],1,sum)
s <- sum(x[ind])
y <- x[-ind]
r <- m[-ind,-ind]
rbind(cbind(r,y),c(y,s))
}
