# Global mutation operator for EA.
#
# @param coords [\code{matrix}]\cr
#   Numeric matrix of city coordinates,
#   rows denote cities.
# @param mutOp [\code{numeric(1)}]\cr
#   Mutation probability from [0,1].
#
# @return [\code{matrix}]
#   Numeric matrix of globally mutated city coordinates.
uniform_mutation = function(coords, mut_op) {
  cities_to_mutate = which(runif(nrow(coords)) < mut_op)
  coords[cities_to_mutate,] = matrix(runif(2 * length(cities_to_mutate)), ncol = 2)
  coords
}

# Local mutation operator of EA.
# Global mutation operator for EA.
#
# @param coords [\code{matrix}]\cr
#   Numeric matrix of city coordinates,
#   rows denote cities.
# @param mut_op [\code{numeric(1)}]\cr
#   Mutation probability from [0,1].
# @param sigma [\code{numeric(1)}]\cr
#   Standard deviation of normal noise.
# @return [\code{matrix}]
#   Numeric matrix of globally mutated city coordinates.
normal_mutation = function(coords, mut_op = 0.1, sigma = 0.0025) {
  cities_to_mutate = which(runif(nrow(coords)) < mut_op)
  ## pmin(pmax(...)) used to ensure we stay in bounds:
  if (length(cities_to_mutate) > 0) {
    delta = matrix(rnorm(2 * length(cities_to_mutate), sd = sigma), ncol = 2)
    coords[cities_to_mutate, ] = pmin(pmax(coords[cities_to_mutate, ] + delta, 0), 1)
  }
  coords
}

# Mating pool generation.
#
# @param pool_size [\code{numeric(1)}]\cr
#   Number of instances to be included to the mating pool.
# @param population [\code{list}]\cr
#   Source population for mating pool generation.
# @param fitness [\code{numeric}]\cr
#   Fitness of members of the population
# @return [\code{list}]
#   List of parents.
create_mating_pool = function(pool_size, population, fitness) {
  mating_pool = list()
  # save population size
  n = length(population)
  for (i in 1:pool_size) {
    # select two random instances from population
    idx = sample(1:n, 2)
    # put member with better fitness value into the mating pool
    mating_pool[[i]] = if (fitness[idx[1]] >= fitness[idx[2]]) {
      population[[idx[1]]]
    } else {
      population[[idx[2]]]
    }
  }
  mating_pool
}

# Round instance (points are placed in the center of the grids).
#
# @param n [\code{numeric(1)}]\cr
#   Number of cells desired (i.e. grid resolution).
# @return [\code{matrix}]
#   Numeric matrix of rounded city coordinates.
round_grid = function(coords, n = 100){
  gr = seq(0, 1, 1/n)
  rnd_grid_pt = apply(coords, 2, function(x) {
    sapply(x, function(y) {
      gr[which.min((y - gr)[(y - gr) >= 0])]
    })
  })

   # avoid outliers outside boundary
  helper = function(x){
    if(all(x!=1)) { y = x + 1 / (2 * n) }
    if(all(x==1)){ y = x - 1 / (2 * n) }
    if((x[2]==1) & (x[1] != 1)) {
      y = c(x[1] + 1 / (2 * n), x[2] - 1 / (2 * n))
    }
    if((x[1]==1) & (x[2] != 1)) {
      y = c(x[1] - 1 / (2 * n), x[2] + 1 / (2 * n))
    }
    return(y)
  }
  t(apply(rnd_grid_pt, 1, helper))
}



#' TSP generating EA.
#'
#' @param fitness_function [\code{function(x, ...)}]\cr
#'   Fitness function used to judge the fitness of a TSP instance.
#'   \code{x} is a numeric matrix with 2 columns, containing
#'   the coordinates of a TSP instance.
#' @param pop_size [\code{integer(1)}]\cr
#'   Number of TSP instances maintained in each population.
#'   Default is 30.
#' @param inst_size [\code{integer(1)}]\cr
#'   Number of cities of each TSP instance.
#'   Default is 50.
#' @param generations [\code{integer(1)}]\cr
#'   Number of generations.
#'   Default is 100L.
#' @param time_limit [\code{integer(1)}]\cr
#'   Time limit in seconds.
#'   Default is 30.
#' @param uniform_mutation_rate [\code{numeric(1)}]\cr
#'   Mutation probability in uniform mutation (in [0,1]).
#' @param normal_mutation_rate [\code{numeric(1)}]\cr
#'   Mutation probability in normal mutation (in [0,1])
#' @param normal_mutation_sd [\code{numeric(1)}]\cr
#'   Standard deviation of normal noise in normal mutation
#' @param cells_round [\code{numeric(1)}]\cr
#'   Grid resolution for rounding
#'   Default is 100.
#' @param rnd [\code{logical(1)}]\cr
#'   Round the coordinates before normal mutation.
#'   Default is \code{TRUE}.
#' @param ... [any]\cr
#'   Not used.
#' @return [\code{list}]
#'   List containing best individual form the last population, its
#'   fitness value, the genrational fitness and the last population.
#'   Default is 50.
#' @export
tsp_generation_ea = function(fitness_function, pop_size = 30L, inst_size = 50L,
  generations = 100L, time_limit = 30L, uniform_mutation_rate,
  normal_mutation_rate, normal_mutation_sd, cells_round = 100L,
  rnd = TRUE, ...) {

  assertFunction(fitness_function, args = c("x"))
  pop_size = convertInteger(pop_size)
  assertInteger(pop_size, len = 1L, any.missing = FALSE, lower = 2L)
  inst_size = convertInteger(inst_size)
  assertInteger(inst_size, len = 1L, any.missing = FALSE, lower = 2L)
  generations = convertInteger(generations)
  assertInteger(generations, len = 1L, any.missing = FALSE, lower = 1L)
  time_limit = convertInteger(time_limit)
  assertInteger(time_limit, len = 1L, any.missing = FALSE, lower = 30)
  assertNumber(uniform_mutation_rate, na.ok = FALSE)
  assertNumber(normal_mutation_rate, na.ok = FALSE)
  assertNumber(normal_mutation_sd, na.ok = FALSE)
  cells_round = convertInteger(cells_round)
  assertInteger(cells_round, len = 1L, any.missing = FALSE)
  assertFlag(rnd, na.ok = FALSE)

  # size of mating pool is half of the population size
  pool_size = round(pop_size / 2)

  ## define overall best problem instance according to fitness value
  overall_best = NULL
  overall_best_fitness = Inf

  ## build initial population randomly by selecting 2 * instSize randoms
  ## (i.e. instSize coordinates) from a R[0,1] distribution
  coords = matrix(runif(pop_size * 2 * inst_size), ncol = 2)
  population = list()
  for (i in 1:pop_size) {
    population_scale = rescale_coords(coords[((i - 1) * inst_size + 1):(i * inst_size), ])
    population[[i]]  = round_grid(population_scale, cells_round)
    if (rnd) {
        population[[i]] = normal_mutation(population[[i]], normal_mutation_rate, normal_mutation_sd)
    }
  }

  start_time = proc.time()[[3]]
  # do the evolutian baby!
  generational_fitness = numeric(generations)
  for (g in 1:generations) {

    # compute fitness value for all instances in current population
    fitness = sapply(population, fitness_function)
    current_time = proc.time()[[3]]
    if (current_time - start_time > time_limit) {
      warning("Time limit reached.")
      break
    }

    current_best = population[[which.min(fitness)]]
    current_best_fitness = min(fitness)
    generational_fitness[g] = current_best_fitness
    if (current_best_fitness < overall_best_fitness) {
      overall_best_fitness = current_best_fitness
      overall_best = current_best
    }
    mating_pool = create_mating_pool(pool_size, population, -fitness)
    message(sprintf("%4i %5.3f %5.3f %3i",
                    g, overall_best_fitness, current_best_fitness,
                    which.max(fitness)))

    # Inspired by the LION 2010 paper we use 1-elitism, i.e. the
    # "best" instance of current population survives with probability 1.
    next_population = vector(length(population), mode = "list")
    next_population[[1]] = current_best
    # 2-tournament selection
    for (k in 2:pop_size) {
      # choose two parents randomly from mating pool
      idx = sample(1:pool_size, 2)
      parent1 = mating_pool[[idx[1]]]
      parent2 = mating_pool[[idx[2]]]

      # build offspring
      offspring = matrix(NA, ncol = 2, nrow = inst_size)
      idx = runif(inst_size) < 0.5
      offspring[idx, ]  = parent1[idx, ]
      offspring[!idx, ] = parent2[!idx, ]

      # mutate
      offspring = uniform_mutation(offspring, uniform_mutation_rate)

      if(rnd) {
        offspring = rescale_coords(offspring)
        offspring = round_grid(offspring, cells_round)
        offspring = normal_mutation(offspring, normal_mutation_rate, normal_mutation_sd)
      } else {
        offspring = normal_mutation(offspring, normal_mutation_rate, normal_mutation_sd)
        offspring = rescale_coords(offspring)
        offspring = round_grid(offspring, cells_round)
      }
      next_population[[k]] = offspring
    }
    ## replace population
    population = next_population
  }
  list(par = current_best, value = current_best_fitness,
    fitness = generational_fitness, last_population = population)
}
