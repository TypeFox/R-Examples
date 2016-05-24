# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( fields)
options( echo=FALSE)
set.seed( 234)
test.for.zero.flag<- 1
n <- 5
m <- 4
mat <- array(rnorm(n*m),c(n,m))
mat2 <- array(rnorm(n*m),c(m,n))
vec <- rnorm(n)
vec2 <- rnorm(n)

test.for.zero( mat2 %*% mat, mat2%d*%mat, tol=1e-8 )

test.for.zero( (diag(vec)%*% mat), (vec%d*%mat), tol=1e-8 )


test.for.zero( diag(vec)%*% vec2, vec%d*%vec2,tol=1e-8)
cat("All done with testing diag multiply", fill=TRUE)
options(echo=TRUE)
