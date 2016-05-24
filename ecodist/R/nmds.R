nmds <- function(dmat, mindim = 1, maxdim = 2, nits = 10, iconf = 0, epsilon = 1e-12, maxit = 500, trace=FALSE, stresscalc="default")

{
# Non-metric multidimensional scaling function
# using the majorization algorithm from
# Borg & Groenen 1997, Modern Multidimensional Scaling.
#
# Sarah Goslee
# 20 December 1999
#
# dmat is a lower-triangular distance matrix.
# mindim is the minimum number of dimensions to calculate.
# maxdim is the maximum number of dimensions to calculate.
# nits is the number of repetitions.
# iconf is the initial configuration to use.
# If iconf is not specified, then a random configuration is used.
# epsilon and maxit specify stopping points.
# Returns a list of configurations (conf)
# and a vector of final stress values (stress),
# along with the cumulative and incremental r^2 for each axis.
# The first nits elements are for the lowest number of dimensions.
#
# stresscalc can be "default" for compatibility with previous versions 
# of ecodist, or "kruskal" to follow Kruskal 1994 method exactly.


nmdscalc <- function(dmat, ndim, iconf, epsilon, maxit, trace)

{

sstress <- function(dmat, cmat, stresscalc)
{
# Calculates the stress-1 function for the original and
# new NMDS configurations (Kruskal 1964).

    if(stresscalc == "default") {
        sstresscalc <- (dmat - cmat) ^ 2
        sstresscalc <- sum(sstresscalc) / sum(cmat ^ 2)

        sqrt(sstresscalc)
    }
    else if(stresscalc == "kruskal") {
        sstresscalc <- (dmat - cmat) ^ 2
        sstresscalc <- sum(sstresscalc) / sum(dmat ^ 2)

        sqrt(sstresscalc)
    }
    else stop("stresscalc has an unknown value\n")

}

# This is the optimization routine for NMDS ordination.
# Use front-end nmds.

n <- (1 + sqrt(1 + 8 * length(dmat))) / 2

if(!is.matrix(iconf)) {
   cat("Using random start configuration \n")
   iconf <- matrix(runif(n * ndim), nrow = n, ncol = ndim)
}

if(dim(iconf)[[2]] != ndim) {
   cat("iconf wrong size: using random start configuration \n")
   iconf <- matrix(runif(n * ndim), nrow = n, ncol = ndim)
}

k <- 0
conf <- iconf
stress2 <- sstress(dmat, dist(iconf), stresscalc)
stress1 <- stress2 + 1 + epsilon

while(k < maxit && abs(stress1 - stress2) > epsilon) {
   stress1 <- stress2

   dmat.full <- full(dmat)
   confd.full <- full(dist(conf))

   confd2.full <- confd.full
   confd2.full[confd.full == 0] <- 1

   b <- dmat.full / confd2.full
   b[confd.full == 0] <- 0
   bsum <- apply(b, 2, sum)
   b <- -1 * b
   diag(b) <- bsum

   conf <- (1 / n) * b %*% conf

   stress2 <- sstress(dmat, dist(conf), stresscalc)

   if(trace) cat(k, ",\t", stress1, "\n")

   k <- k + 1
}

list(conf = conf, stress = stress1)

}


conf <- list(1:((maxdim - mindim + 1) * nits))
stress <- list(1:((maxdim - mindim + 1) * nits))
r2 <- list(1:((maxdim - mindim + 1) * nits))

k <- 1

for(i in mindim:maxdim) {
   if(trace) cat("Number of dimensions: ", i, "\n")
   for(j in 1:nits) {
      if(trace) cat("Iteration number: ", j, "\n")
      nmdsr <- nmdscalc(dmat, ndim = i, iconf, epsilon, maxit, trace)
      conf[[k]] <- nmdsr$conf
      stress[[k]] <- nmdsr$stress
      r2[[k]] <- cor(dmat, dist(nmdsr$conf)) ^ 2
      k <- k + 1
   }
}

list(conf = conf, stress = unlist(stress), r2 = unlist(r2))

}

