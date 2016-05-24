# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

permutationTest <- function (data, statistic, R = 9999,
                             alternative = "two.sided",
                             resampleColumns = NULL,
                             args.stat = NULL,
                             seed = NULL,
                             sampler = samp.permute,
                             label = NULL, statisticNames = NULL,
                             block.size = 100, trace = FALSE,
                             tolerance = .Machine$double.eps ^ 0.5)
{
  # Permutation test

  # Args:
  #   data:      vector, matrix, or data frame.
  #              Let 'nObs' = length of vector, or nrow otherwise.
  #   statistic: a function, or expression (e.g. mean(data, trim = .2)
  #              If data is a data frame, can refer to variables in it.
  #              This may be a vector; let 'd' be its length.
  #   R:         number of replications
  #   resampleColumns: column indices or names; if specified, then only
  #              these columns are permuted (data must be a matrix or
  #              data frame).
  #   alternative: one of "two.sided", "greater", or "less"
  #   args.stat: additional arguments to pass to the function
  #   seed:      old value of .Random.seed, or argument to set.seed
  #   sampler:   a function for resampling, see help(samp.permute)
  #   label:     used for labeling plots
  #   statisticNames: names used for printing, character vector of length 'd'
  #   block.size: replicates are done 'block.size' at a time
  #   trace:     logical, if TRUE an indication of progress is printed.
  #   tolerance: numerical tolerance when computing P-values; smaller
  #              differences between replicated & observed are considered equal.

  Call <- match.call()
  resampleFun <-
    .resampleMakeFunction(data, statistic,
                          substitute(data), substitute(statistic), args.stat,
                          resampleColumns = resampleColumns)
  result <- resample(data, resampleFun, sampler = sampler, R = R,
                     seed = seed,
                     statisticNames = statisticNames,
                     block.size = block.size, trace = trace, call = Call)
  result$stats <- .PermutationStats(result, alternative = alternative,
                                    tolerance = tolerance)
  class(result) <- c("permutationTest", "resample")
  result
}

# permutationTest(y, cor(y, x))
# permutationTest(myData, cor(y, x), resampleColumns = "y")


# print.resample should suffice


if(FALSE) {
  x9 <- c(2:3, 1:7)
  xDF <- data.frame(X = x9, Y = 2*x9)
  mystat <- function(x) sum(x * seq_along(x))
  mycor <- function(x) cor(x, seq_len(if(length(dim(x)) == 2) nrow(x) else length(x)))



  source("~/resample/R/permutationTest.R")

  ### statistic by name
  # base case: data by name, statistic is function by name
  permutationTest(x9, mystat, R=99)

  # data expression
  permutationTest((x9), mystat, R=99)

  # inline function
  permutationTest(x9, function(z) mystat(z), R=99)

  # data frame,
  permutationTest(xDF, mycor, R=99)

  # data expression, data frame
  permutationTest((xDF), mycor, R=99)

  # data expression, matrix
  permutationTest(as.matrix(xDF), mycor, R=99)


  ### statistic expression
  # data by name
  permutationTest(x9, mystat(x9), R=99)

  # data as expression, refer to 'data'
  permutationTest((x9), mystat(data), R=99)

  # data frame
  permutationTest(xDF, mystat(X), R=99)

  # data frame expression
  permutationTest((xDF), mystat(X), R=99)

  # See if results reproduce
  temp <- .Last.value
  .Random.seed <- temp$seed
  all.equal(temp, eval(temp$call))

  # resampleColumns
  permutationTest(xDF, cor(X, Y), resampleColumns = "X", R=99)

  source("~/resample/R/permutationTest.R")
}
