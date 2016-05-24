library("aroma.core")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# True CN states
stateFcn <- function(data, ...) {
  x <- data$x
  states <- integer(length(x))
  states[100 <=x & x <= 150] <- -1L
  states[320 <=x & x <= 400] <- +1L
  states
}

stateFcn <- function(data, ...) {
  x <- data$x
  states <- rep("neutral", time=length(x))
  states[100 <=x & x <= 150] <- "loss"
  states[320 <=x & x <= 400] <- "gain"
  states
}

# Number of loci
J <- 500

eps <- rnorm(J, sd=1/2)
mu <- double(J)
x <- 1:J

levels <- c("neutral"=0, "loss"=-1, "gain"=+1);

cn <- RawCopyNumbers(eps, x=x)
cn$state <- stateFcn
cn$cn <- cn$cn + 0.8*levels[cn$state]
print(cn)


cn2 <- extractSubset(cn, subset=xSeq(cn, by=5))
print(cn2)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot along genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, ylim=c(-3,3))
title(main="Complete and subsetted loci")
points(cn2, col="red", pch=176, cex=2)

legend("topright", pch=c(19,176), col=c("#999999", "red"), sprintf(c("raw [n=%d]", "every 5th [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cn2))), bty="n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Binned smoothing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, col="#999999", ylim=c(-3,3))
title(main="Binned smoothing by State")

cnSa <- binnedSmoothingByField(cn, by=3, field="state")
lines(cnSa, col="blue")
points(cnSa, col="blue")

cnSb <- binnedSmoothingByField(cn, by=9, field="state")
lines(cnSb, col="red")
points(cnSb, col="red")

legend("topright", pch=19, col=c("#999999", "blue", "red"), sprintf(c("raw [n=%d]", "Bin(w=3) [n=%d]", "Bin(w=9) [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cnSa), nbrOfLoci(cnSb))), bty="n")
