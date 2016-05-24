# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

bootstrap <- function(data, statistic, R = 10000, args.stat = NULL,
                      seed = NULL,
                      sampler = samp.bootstrap,
                      label = NULL, statisticNames = NULL,
                      block.size = 100, trace = FALSE)
{
  # Nonparametric bootstrap.
  #
  # Args:
  #   data:      vector, matrix, or data frame.
  #              Let 'nObs' = length of vector, or nrow otherwise.
  #   statistic: a function, or expression (e.g. mean(data, trim = .2)
  #              If data is a data frame, can refer to variables in it.
  #              This may be a vector; let 'd' be its length.
  #   R:         number of replications
  #   args.stat: additional arguments to pass to the function
  #   seed:      old value of .Random.seed, or argument to set.seed
  #   sampler:   a function for resampling, see help(samp.bootstrap)
  #   label:     used for labeling plots
  #   statisticNames: names used for printing, character vector of length 'd'
  #   block.size: replicates are done 'block.size' at a time
  #   trace:     logical, if TRUE an indication of progress is printed.

  Call <- match.call()
  resampleFun <-
    .resampleMakeFunction(data, statistic,
                          substitute(data), substitute(statistic), args.stat)
  result <- resample(data, resampleFun, sampler = sampler, R = R,
                     seed = seed,
                     statisticNames = statisticNames,
                     block.size = block.size, trace = trace, call = Call)
  result$stats <- .BootStats(result)
  class(result) <- c("bootstrap", "resample")
  result
}
# TODO: support group

# print.resample should suffice
# print.bootstrap <- function(x, ...) {
#   cat0n("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"))
#   catn("Replications:", x$R)
#   catn("\nSummary Statistics:")
#   print(x$stats, ...)
#   invisible(x)
# }


if(FALSE) {
  x9 <- 1:9
  xDF <- data.frame(X = x9, Y = 2*x9)

  source("~/resample/R/bootstrap.R")

  ### statistic by name
  # base case: data by name, statistic is function by name
  bootstrap(x9, mean, R=100)

  # data expression
  bootstrap((x9), mean, R=100)

  # args.stat
  bootstrap(x9, mean, args.stat = list(trim = .25), R=100)

  # inline function
  bootstrap(x9, function(z) mean(z), R=100)

  # data frame,
  bootstrap(xDF, colMeans, R=100)

  # data expression, data frame
  bootstrap((xDF), colMeans, R=100)

  # data expression, matrix
  bootstrap(as.matrix(xDF), colMeans, R=100)


  ### statistic expression
  # data by name
  bootstrap(x9, mean(x9), R=100)

  # data as expression, refer to 'data'
  bootstrap((x9), mean(data), R=100)

  # data frame
  bootstrap(xDF, mean(X), R=100)

  # data frame expression
  bootstrap((xDF), mean(X), R=100)

  # See if results reproduce
  temp <- .Last.value
  .Random.seed <- temp$seed
  all.equal(temp, eval(temp$call))

  source("~/resample/R/bootstrap.R")
}
