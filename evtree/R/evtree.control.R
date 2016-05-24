evtree.control <- function(minbucket = 7L, minsplit = 20L, maxdepth = 9L,
  niterations = 10000L, ntrees = 100L, alpha = 1,
  operatorprob = list(pmutatemajor = 0.2, pmutateminor = 0.2, pcrossover = 0.2, psplit = 0.2, pprune = 0.2),
  seed = NULL, ...)
{
  minbucket <- as.integer(minbucket)
  if(minbucket < 1L) {
    warning("parameter \"minbucket\" must be positive, default used instead")
    minbucket <- 7L
  }

  minsplit <- as.integer(minsplit)
  if(minsplit < 2 * minbucket) {
    if(!missing(minsplit)) warning("parameter \"minsplit\" must be at least twice as large as \"minbucket\", changed")
    minsplit <- 2 * minbucket
  }

  maxdepth <- as.integer(maxdepth) 
  if(maxdepth > 15) {
    warning("computations may take extremely long for \"maxdepth\" > 15 (or even be infeasible)")
  }else if(maxdepth < 1) {
    warning("parameter \"maxdepth\" must be equal or larger than 1, default used")
    maxdepth <- 9
  }

  if(!is.integer(niterations)) niterations <- as.integer(niterations)
  if(niterations < 100L) {
    warning("computations may be unreliable for \"niterations\" < 100")
  }

  if(!is.integer(ntrees)) ntrees <- as.integer(ntrees)
  if(ntrees < 2L) {
    warning("parameter \"ntrees\" must be equal or larger than 2, default used")	
    ntrees <- 100
  }else if(ntrees < 10L) {
    warning("computations may be unreliable for \"ntrees\" < 10")
  }

  alpha <- as.numeric(alpha)
  if(alpha < 0) {
    warning("parameter \"alpha\" must be equal or larger than 0, default used")
    alpha <- 1
  }

  if(is.null(seed)) seed <- as.integer(runif(1, max = 2^16))
  if(!is.integer(seed)) seed <- as.integer(seed)
  if(seed < -1L) {
    warning("parameter \"seed\" must be non-negative, default used instead")
    seed <- -1L
  }

  op <- c(pmutatemajor = 0.2, pmutateminor = 0.2, pcrossover = 0.2, psplit = 0.2, pprune = 0.2)
  operatorprob <- unlist(operatorprob)
  if(is.null(names(operatorprob)) & length(operatorprob) == 5L) names(operatorprob) <- names(op)
  if(is.null(names(operatorprob))) {
    warning("incorrect specification of \"operatorprob\", default used")
    operatorprob <- op
  }
  op[names(operatorprob)] <- operatorprob
  operatorprob <- op/sum(op)*100

  list(
    minbucket = minbucket,
    minsplit = minsplit,
    maxdepth = maxdepth,
    niterations = niterations,
    ntrees = ntrees,
    alpha = alpha,
    operatorprob = operatorprob,
    seed = seed
  )
}

