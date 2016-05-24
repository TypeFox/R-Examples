## Goodman (1991), Table 17.1 (p. 1097)
# See also ?rc

library(logmult)
data(criminal)

model <- rc(criminal, start=NA)

model$assoc # These are the phi (.07), mu and nu
model$assoc$row[,1,1] * model$assoc$phi[1,1] # These are the mu'
model$assoc$col[,1,1] * model$assoc$phi[1,1] # These are the nu'

stopifnot(round(model$assoc$phi[1,1], d=2) == 0.07)
# Scores have reversed signs compared to the orignal article
stopifnot(isTRUE(all.equal(round(model$assoc$row[,1,1], d=2),
                           c(1.26, 0.82, -0.56, -1.21),
                           check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(round(model$assoc$col[,1,1], d=2),
                           c(-1.44, -1.30, -0.33, 1.00, 0.89),
                           check.attributes=FALSE)))
