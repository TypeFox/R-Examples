# From SamplerCompare, (c) 2010 Madeleine Thompson

# This script checks ensures that some of the example distributions
# have gradient functions that match their log density functions.

library(SamplerCompare)

dists <- list(N2weakcor.dist, N4poscor.dist, N4negcor.dist, schools.dist)

for (dist in dists) {
  check.dist.gradient(dist, dist$mean + 0.1)
}
