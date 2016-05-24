## Goodman (1979), scores not reported but present in LEM examples
# as "GOO91_1G" (where scores reported here come from)
# No symmetric RC scores seem to be reported in a published paper

library(logmult)
data(occupationalStatus)

# We use a simpler (collapsed) table
occupationalStatus[5,]<-colSums(occupationalStatus[5:6,])
occupationalStatus[,5]<-rowSums(occupationalStatus[,5:6])
occupationalStatus <- occupationalStatus[-6,-6]

model <- rc(occupationalStatus, diagonal=TRUE, symmetric=TRUE, weighting="none", start=NA)

stopifnot(round(model$assoc$phi[1,1], d=3) == 6.10)
stopifnot(isTRUE(all.equal(round(c(model$assoc$row), d=3),
                           c(0.532, 0.438, 0.206, -0.031,
                            -0.216, -0.426, -0.503))))

model <- rc(occupationalStatus, diagonal=TRUE, symmetric=TRUE, weighting="uniform", start=NA)

stopifnot(round(model$assoc$phi[1,1], d=3) == 0.871)
stopifnot(isTRUE(all.equal(round(c(model$assoc$row), d=3),
                           c(1.409, 1.159, 0.544, -0.082,
                            -0.571, -1.127, -1.332))))

# Marginal weighted scores from ass2 in LEM are weird: we need
# to run cor(3) on the fitted table to get correct values, which is not
# possible for symmetric scores
