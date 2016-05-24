# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

jackknife <- function(data, statistic, args.stat = NULL,
                      label = NULL, statisticNames = NULL,
                      trace = FALSE)
{
  # Basic jackknife
  #
  # Args:
  #   data:      vector, matrix, or data frame.
  #              Let 'nObs' = length of vector, or nrow otherwise.
  #   statistic: a function, or expression (e.g. mean(data, trim = .2)
  #              If data is a data frame, can refer to variables in it.
  #              This may be a vector; let 'd' be its length.
  #   args.stat: additional arguments to pass to the function
  #   label:     used for labeling plots
  #   statisticNames: names used for printing, character vector of length 'd'
  #   block.size: replicates are done 'block.size' at a time
  #   trace:     logical, if TRUE an indication of progress is printed.

  Call <- match.call()
  resampleFun <-
    .resampleMakeFunction(data, statistic,
                          substitute(data), substitute(statistic), args.stat)
  dimData <- dim(data)
  n <- IfElse(is.null(dimData), length(data), dimData[1])

  result <- resample(data, resampleFun,
                     sampler = function(n, R) -matrix(1:n, 1), R = n,
                     statisticNames = statisticNames,
                     block.size = n, trace = trace,
                     call = match.call())
  result$stats <- .JackknifeStats(result)
  class(result) <- c("bootstrap", "resample")
  result
}



if(FALSE) {
  x9 <- (1:9)^2
  xDF <- data.frame(X = x9, Y = 2*x9)

  source("~/resample/R/jackknife.R")

  ### statistic by name
  # base case: data by name, statistic is function by name
  jackknife(x9, mean)

  # data expression
  jackknife((x9), mean)

  # args.stat
  jackknife(x9, mean, args.stat = list(trim = .25))

  # inline function
  jackknife(x9, function(z) mean(z))

  # data frame,
  jackknife(xDF, colMeans)

  # data expression, data frame
  jackknife((xDF), colMeans)

  # data expression, matrix
  jackknife(as.matrix(xDF), colMeans)


  ### statistic expression
  # data by name
  jackknife(x9, mean(x9))

  # data as expression, refer to 'data'
  jackknife((x9), mean(data))

  # data frame
  jackknife(xDF, mean(X))

  # data frame expression
  jackknife((xDF), mean(X))

  # See if results reproduce
  temp <- .Last.value
  .Random.seed <- temp$seed
  all.equal(temp, eval(temp$call))

  source("~/resample/R/jackknife.R")
}
