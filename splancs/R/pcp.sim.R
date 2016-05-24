# Copyright Giovanni Petris <GPetris@uark.edu> and Nicolas Picard
# <nicolas.picard@cirad.fr> 2001-2003
#
###
### Generate a Poisson Cluster Process
###
pcp.sim <- function(rho, m, s2, region.poly, larger.region=NULL,
 vectorise.loop=TRUE) {
  ## rho: intensity of the parent process
  ## m: average number of offsprings per parent
  ## s2: variance of location of offsprings relative to
  ##   their parent
  ## region.poly: a polygon defining the region in which
  ##   the process is to be generated
  ## larger.region: a rectangle containing the region of interest
  ##   given in the form (xl,xu,yl,yu)
  ## vectorise.loop: if TRUE use vectorised code by Nicolas Picard
  ##   otherwise use original loop and rbind
  if (is.null(larger.region))
    larger.region <- as.vector(apply(sbox(region.poly), 2, range))
  sim.events <- c(0,0)
  ## 1. Generate the parents on [xl,xu]x[yl,yu]
  n <- rpois(1,lambda=rho*(larger.region[2]-larger.region[1])*
             (larger.region[4]-larger.region[3]))
  parents <- cbind(runif(n,larger.region[1],larger.region[2]),
                   runif(n,larger.region[3],larger.region[4]))
  ## 2. Generate the children
  sd <- sqrt(s2)
  if (!vectorise.loop) {
    for (j in 1:n) {
      num.child <- rpois(1,lambda=m)
      for (k in 1:num.child) {
        new.child <- parents[j,]+rnorm(2,0,sd=sd)
        sim.events <- rbind(sim.events,new.child)
      }
    }
    sim.events <- sim.events[-1,]  
  } else {
# vectorisation contributed by Nicolas Picard
    num.child <- rpois(n, lambda = m)
    num.tot <- sum(num.child)
    sim.events <- matrix(c(rnorm(num.tot, 0, sd) + rep(parents[, 1],
      num.child), rnorm(num.tot, 0, sd) + rep(parents[ ,2],
      num.child)), num.tot)
  }
  ## return only the events within the region of interest
  pip(as.points(sim.events),region.poly)
}

