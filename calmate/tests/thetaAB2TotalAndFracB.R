library("calmate")

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate")
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path)

data <- thetaAB2TotalAndFracB(theta)
str(data)

theta2 <- totalAndFracB2ThetaAB(data)
str(theta2)

stopifnot(all.equal(theta2, theta))

