# program spuRs/resources/scripts/ev_web.r

popsize <- 16 # size of population, must be even
r_mutate <- 0.25 # proportion of each generation mutated
# set up plotting
opar <- par(mfrow=c(4,4), mar=c(0,0,0,0), oma=c(0,0,0,0))
# set up initial population of webs, each consisting of a
# single strand
population <- list()
for (i in 1:popsize) {
  population <- c(population, list(runif(2, 0, 2*pi)))
}
# main loop
continue <- TRUE # flag to continue searching
n_gen <- 1 # num of generations in next batch
total_gen <- 0 # total generations so far
while (continue) {
  for (gen in 1:n_gen) {
    # mutate
    for (i in 1:popsize) {
      if (runif(1) < r_mutate) {
        population[[i]] <- mutate(population[[i]])
      }
    }
    # recombine
    x <- sample(1:popsize, popsize/2)
    y <- (1:popsize)[-x]
    for (i in 1:(popsize/2)) {
      population <- c(population, recombine(population[[x[i]]],
                                            population[[y[i]]]))
    }
    # select
    fitness <- rep(0, 2*popsize)
    for (i in 1:(2*popsize)) {
      fitness[i] <- efficiency(population[[i]])
    }
    fittest <- order(fitness, decreasing=TRUE)
    population <- population[fittest[1:popsize]]
    cat(".") # so we know something is happening
  }
  # report
  total_gen <- total_gen + n_gen
  cat("\nbest efficiency after", total_gen,
      "generations is", fitness[fittest[1]], "\n")
  for (i in 1:popsize) webdraw(population[[i]])
  # get num of generations in next batch, zero to stop
  n_gen <- as.numeric(readline("number of extra generations: "))
  continue <- n_gen > 0
}
par <- opar
