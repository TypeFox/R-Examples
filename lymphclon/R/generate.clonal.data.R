#!/usr/bin/Rscript

generate.clonal.data <- function(
n = 2e7, # typical number of B cell clones in an human -- typically 20 million
num.cells.taken.vector = c(2e3, 5e3, 1e4, 2e4, 5e4, 5e4),
read.count.per.replicate.vector = rep(2e4, length(num.cells.taken.vector)),
clonal.distribution.power = -sqrt(2),
pcr.noise.type = 'pareto',
pcr.pareto.location = 1,
pcr.pareto.shape = 1,
pcr.lognormal.meanlog = 0,
pcr.lognormal.sdlog = 1)
{
pcr.rpareto.func <- function(n, location, shape) {
  # taken from VGAM 0.92
  ans <- location / runif(n)^(1/shape)
  ans[location <= 0] <- NaN
  ans[shape    <= 0] <- NaN
  ans
}
unnorm.clone.prob <- (1:n) ^ clonal.distribution.power
true.clone.prob <- unnorm.clone.prob / sum(unnorm.clone.prob)
true.clonality <- sum(true.clone.prob * true.clone.prob)
num.replicates <- length(num.cells.taken.vector)
replicates <- matrix(0, n, num.replicates)
replicate.errs <- matrix(0, n, num.replicates)
replicate.squared.errs <- rep(NA, num.replicates)

if (pcr.noise.type == 'pareto') {
  get.readcount.given.cellcount <- function(x) {
    ifelse(x > 0, sum(abs(pcr.rpareto.func(n = x, location = pcr.pareto.location, shape = pcr.pareto.shape))), 0)
  }
} else if (pcr.noise.type == 'lognormal') { # lognormal
  get.readcount.given.cellcount <- function(x) {
    ifelse(x > 0, sum(abs(rlnorm(n = x, meanlog = pcr.lognormal.meanlog, sdlog = pcr.lognormal.sdlog))), 0)
  }
} else {
  print(sprintf('Unknown pcr noise type: %s', pcr.noise.type))
}

for (i in 1:num.replicates)
{

  read.count.per.replicate <- read.count.per.replicate.vector[i]
  num.cells.taken <- num.cells.taken.vector[i]

  num.obs.clones.frac <- num.cells.taken / n
  sample.of.cells <- rmultinom(n = 1, size = num.cells.taken, prob = true.clone.prob)
  raw.sample.of.reads <- sapply(X = sample.of.cells, FUN = get.readcount.given.cellcount)
  sample.of.reads <- rpois(
    n = length(sample.of.cells),
    lambda = raw.sample.of.reads * (read.count.per.replicate / sum(raw.sample.of.reads)))
  # table(sample.of.reads)
  # return the counts, rather than the abundance proportions
  replicates[, i] <- sample.of.reads # / sum(sample.of.reads)
  replicate.errs[, i] <- sample.of.reads / sum(sample.of.reads) - true.clone.prob
  replicate.squared.errs[i] <- sum(replicate.errs[, i] ^ 2)
  grahm.errs <- t(replicate.errs) %*% replicate.errs
}

simulated.data <- 
  list(read.count.matrix = replicates, 
       true.clone.prob = true.clone.prob, 
       true.clonality = true.clonality,
       replicate.squared.errs = replicate.squared.errs,
       grahm.errs = grahm.errs)

return(simulated.data)
}
