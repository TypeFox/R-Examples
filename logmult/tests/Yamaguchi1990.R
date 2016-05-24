## Yamaguchi (1990), Table 5, p. 202, and Table 6B, p. 205

library(logmult)
data(ocg1973)

# Simple symmetric RC(1) model ("Null skew-symmetry")
rc.model <- rc(ocg1973, diagonal=TRUE, symmetric=TRUE, weighting="none")

# Reported phi is slightly different, coefficients agree
rc.model
stopifnot(isTRUE(all.equal(round(rc.model$assoc$row[,1,1], d=2),
                           c(0.69, 0.3, -0.1, -0.33, -0.56),
                           check.attributes=FALSE)))


# Note model does not always converge, several attempts may be needed
# Here we set known starting values to be sure it works
set.seed(5)
model <- yrcskew(ocg1973, nd.symm=1, nd.skew=1, diagonal=TRUE, weighting="none")

# We do not get the same results as the author, but the smaller deviance
# indicates a better fit in our version (!)
model


summary(model)

# These scores cannot be checked since the author reports a higher deviance:
# they are simply here to ensure no regressions appear
stopifnot(isTRUE(all.equal(c(round(model$assoc.yrcskew$row, d=3)),
                           c(-0.421, -0.332, -0.166, 0.097, 0.822))))
stopifnot(isTRUE(all.equal(c(round(model$assoc$row, d=3)),
                             c(0.851, 0.049, -0.246, -0.331, -0.322))))
