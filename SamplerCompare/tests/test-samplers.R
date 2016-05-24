# From SamplerCompare, (c) 2010 Madeleine Thompson

# This script makes sure that all the samplers exported in NAMESPACE
# can draw a 500-state sample from a simple 2D Gaussian implemented
# in both R and C, and that they generate results that are not horribly
# wrong.

# Since this has to run every time the package is checked, the chain
# length must be short.

library(SamplerCompare)

all.samplers <- list(
  shrinking.rank.sample, nonadaptive.crumb.sample,
  stepout.slice.sample, interval.slice.sample, hyperrectangle.sample,
  nograd.hyperrectangle.sample, arms.sample, cov.match.sample,
  multivariate.metropolis.sample, univar.metropolis.sample,
  adaptive.metropolis.sample, univar.eigen.sample, cheat.univar.eigen.sample,
  oblique.hyperrect.sample, cheat.oblique.hyperrect.sample)

# Use a second core if synchronicity package is available.  In production,
# we would only do this if we were on a multi-core system, but it's
# better to lose some efficiency here and guarantee the code path is
# tested.  If we wanted to know the actual number of cores, the best
# way to find it is parallel:::detectCores().

has.synchronicity <- 'synchronicity' %in% installed.packages()[,'Package']
if (has.synchronicity) {
  cores <- 2
} else {
  cores <- 1
}

RS <- compare.samplers(700, list(N2weakcor.dist), all.samplers,
                       burn.in=0.7, cores=cores)
stopifnot(max(RS$err)<1)

# A version of N2weakcor.dist implemented in C.

N2weakcor.dist.C <- make.c.dist(2, 'Gauss2-C', 'Gauss2_log_dens',
  c(N2weakcor.dist$mean, 0.8), mean=N2weakcor.dist$mean, cov=N2weakcor.dist$cov)

RS <- compare.samplers(700, list(N2weakcor.dist.C), all.samplers,
                       burn.in=0.7, cores=cores)
stopifnot(max(RS$err)<1)
