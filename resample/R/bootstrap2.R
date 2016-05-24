# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

bootstrap2 <- function(data, statistic,
                       treatment, data2 = NULL,
                       R = 10000, ratio = FALSE, args.stat = NULL,
                       seed = NULL,
                       sampler = samp.bootstrap,
                       label = NULL, statisticNames = NULL,
                       block.size = 100, trace = FALSE)
{
  # Two-sample nonparametric bootstrap.
  # Specify either (data, data2) or (data, treatment).
  #
  # Args:
  #   data:      vector, matrix, or data frame.
  #              Let 'nObs' = length of vector, or nrow otherwise.
  #   statistic: a function, or expression (e.g. mean(data, trim = .2)
  #              If data is a data frame, can refer to variables in it.
  #              This may be a vector; let 'd' be its length.
  #   treatment: a vector of length nObs, with two unique values.
  #   data2:     like data.
  #   R:         number of replications.
  #   ratio:     logical, if FALSE then bootstrap the difference between
  #              data and data2 (or first and second treatment).
  #              If TRUE then bootstrap the ratio.
  #   args.stat: additional arguments to pass to the function.
  #   seed:      old value of .Random.seed, or argument to set.seed.
  #   sampler:   a function for resampling, see help(samp.bootstrap).
  #   label:     used for labeling plots.
  #   statisticNames: names used for printing, character vector of length 'd'.
  #   block.size: replicates are done 'block.size' at a time.
  #   trace:     logical, if TRUE an indication of progress is printed.

  Call <- match.call()

  dimData <- dim(data)
  n <- IfElse(is.null(dimData), length(data), dimData[1])
  stopifnot(length(dimData) <= 2)

  resultsBoth <- vector("list", 2)
  if(is.null(data2)) {                  # use treatment
    if(is.data.frame(data))
      treatment <- eval(substitute(treatment),
                        envir = data, enclos = parent.frame())
    treatmentInds <- split(1:n, treatment)
    treatmentNames <- names(treatmentInds)
    if(length(treatmentInds) != 2)
      stop("treatment must have 2 unique values, observed: ",
           paste(treatmentNames, collapse = " "))
    for(k in 1:2) {
      # Create a smaller data set corresponding to this treatment
      dataK <- IfElse(length(dimData) == 2,
                      data[treatmentInds[[k]], , drop = FALSE],
                      data[treatmentInds[[k]]])
      resampleFun <-
        .resampleMakeFunction(dataK, statistic,
                              substitute(data), substitute(statistic),
                              args.stat)
      resultsBoth[[k]] <-
        resample(dataK, resampleFun, sampler = sampler, R = R,
                 seed = IfElse(k == 1, seed, .Random.seed),
                 statisticNames = statisticNames,
                 block.size = block.size, trace = trace, call = Call)
    }
  } else {  # data and data2
    if(length(dim(data2)) != length(dimData) ||
       class(data2) != class(data))
      stop("data and data2 must have similar structure")
    treatmentNames <- c(IfElse(is.name(substitute(data)),
                               substitute(data), "data"),
                        IfElse(is.name(substitute(data2)),
                               substitute(data2), "data2"))
    for(k in 1:2) {
      resampleFun <-
        .resampleMakeFunction(IfElse(k == 1, data, data2), statistic,
                              as.name("data"), substitute(statistic), args.stat)
      resultsBoth[[k]] <-
        resample(IfElse(k == 1, data, data2),
                 resampleFun, sampler = sampler, R = R,
                 seed = IfElse(k == 1, seed, .Random.seed),
                 statisticNames = statisticNames,
                 block.size = block.size, trace = trace, call = Call)
    }
  }
  for(k in 1:2) {
    resultsBoth[[k]]$stats <- .BootStats(resultsBoth[[k]])
    resultsBoth[[k]]$call <- paste("bootstrap for data set", k)
    class(resultsBoth[[k]]) <- c("bootstrap", "resample")
  }
  names(resultsBoth) <- treatmentNames
  if(!.CheckCompatibleObserved(resultsBoth))
    return(resultsBoth)

  # combine results
  op <- get(IfElse(ratio, "/", "-"))
  result <- list(observed =
                 op(resultsBoth[[1]]$observed, resultsBoth[[2]]$observed),
                 replicates =
                 op(resultsBoth[[1]]$replicates, resultsBoth[[2]]$replicates),
                 n = c(resultsBoth[[1]]$n, resultsBoth[[2]]$n),
                 p = resultsBoth[[1]]$p,
                 R = R,
                 seed = resultsBoth[[1]]$seed,
                 call = Call,
                 resultsBoth = resultsBoth,
                 ratio = ratio)
  result$failures <- unique(c(resultsBoth[[1]]$failures,
                              resultsBoth[[2]]$failures))
  statisticNames <- .StatisticNames2(statisticNames, treatmentNames, ratio,
                                     resultsBoth[[1]]$observed, Call)
  names(result$observed) <- statisticNames
  colnames(result$replicates) <- statisticNames
  result$stats <- .BootStats(result)
  class(result) <- c("bootstrap2", "bootstrap", "resample")
  result
}

# Note difference:
#   bootstrap(x, mean(x))
#   bootstrap2(x1, data2 = x2, mean(data))

.CheckCompatibleObserved <- function(x){
  # x is a pair of objects from bootstrap2 or permutationTest2.
  # Check that observed statistics have the same names
  # return(TRUE) if they do.
  if(length(x[[1]]$observed) == length(x[[2]]$observed) &&
     identical(names(x[[1]]$observed), names(x[[2]]$observed)))
    return(TRUE)
  warning("Statistics returned from the two samples are not compatible,",
          "will return results so far, without combining them.")
  print(all.equal(x[[1]]$observed, x[[2]]$observed))
  return(FALSE)
}

.StatisticNames2 <- function(statisticNames, treatmentNames, ratio,
                             observed1, call){
  # Create statisticNames like "mean: data1 - data2"
  # statisticNames correspond to one sample
  if(is.null(statisticNames))
    statisticNames <- names(observed1)
  p <- length(observed1)
  if(is.null(statisticNames) && is.call(call) && is.name(call$statistic))
    statisticNames <- IfElse(p == 1, as.character(call$statistic),
                             paste0(as.character(call$statistic), 1:p))
  if(is.null(statisticNames))
    statisticNames <- paste0("stat", 1:p)
  temp <- which(statisticNames == "")
  if(length(temp))
    statisticNames[temp] <- paste0("stat", temp)
  statisticNames <- paste0(statisticNames, ": ",
                           paste(treatmentNames,
                                 collapse = IfElse(ratio, "/", "-")))
}

# print.bootstrap2 <- function(x, ...) {
#   cat0n("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"))
#   catn("Replications:", x$R)
#   catn("Two groups, sample sizes are", x$n[1], x$n[2])
#   catn("\nSummary Statistics for the", IfElse(x$ratio, "ratio", "difference"),
#        "between samples 1 and 2:")
#   print(x$stats, ...)
#   invisible(x)
# }
# same as print.permutatonTest2


if(FALSE) {
  x9 <- 1:9
  xDF <- data.frame(X = x9, Y = 2*x9)
  t9 <- letters[c(1, 1, 1, 2, 1, 2, 1, 2, 1)]
  xDF1 <- xDF[t9 == "a", ]
  xDF2 <- xDF[t9 == "b", ]
  x91 <- x9[t9 == "a"]
  x92 <- x9[t9 == "b"]

  source("~/resample/R/bootstrap2.R")

  ##### treatment
  ### statistic by name
  # base case: data by name, statistic is function by name
  temp <- bootstrap2(treatment = t9, x9, mean, R=100)
  temp

  sapply(temp$resultsBoth, class)
  temp$resultsBoth

  # data expression
  bootstrap2(treatment = t9, (x9), mean, R=100)

  # args.stat
  bootstrap2(treatment = t9, x9, mean, args.stat = list(trim = .25), R=100)

  # inline function
  bootstrap2(treatment = t9, x9, function(z) mean(z), R=100)

  # data frame
  bootstrap2(treatment = t9, xDF, colMeans, R=100)

  # data expression, data frame
  bootstrap2(treatment = t9, (xDF), colMeans, R=100)

  # data expression, matrix
  bootstrap2(treatment = t9, as.matrix(xDF), colMeans, R=100)


  ### statistic expression
  # data by name
  bootstrap2(treatment = t9, x9, mean(x9), R=100)

  # data as expression, refer to 'data'
  bootstrap2(treatment = t9, (x9), mean(data), R=100)

  # data frame
  bootstrap2(treatment = t9, xDF, mean(X), R=100)

  # data frame expression
  bootstrap2(treatment = t9, (xDF), mean(X), R=100)

  # See if results reproduce
  temp <- .Last.value
  .Random.seed <- temp$seed
  all.equal(temp, eval(temp$call))


  ##### data and data2
  ### statistic by name
  # base case: data by name, statistic is function by name
  temp <- bootstrap2(x91, data2 = x92, mean, R=100)
  temp

  sapply(temp$resultsBoth, class)
  temp$resultsBoth

  # data expression
  bootstrap2((x91), data2 = (x92), mean, R=100)

  # args.stat
  bootstrap2(x91, data2 = x92, mean, args.stat = list(trim = .25), R=100)

  # inline function
  bootstrap2(x91, data2 = x92, function(z) mean(z), R=100)

  # data frame
  bootstrap2(xDF1, data2 = xDF2, colMeans, R=100)

  # data expression, data frame
  bootstrap2((xDF1), data2 = (xDF2), colMeans, R=100)

  # data expression, matrix
  bootstrap2(as.matrix(xDF1), data2 = as.matrix(xDF2), colMeans, R=100)


  ### statistic expression
  # data by name
  bootstrap2(x91, data2 = x92, mean(data))  # use 'data' in this case

  # data as expression, refer to 'data'
  bootstrap2((x91), data2 = (x92), mean(data), R=100)

  # data frame
  bootstrap2(xDF1, data2 = xDF2, mean(X), R=100)

  # data frame expression
  bootstrap2((xDF1), data2 = (xDF2), mean(X), R=100)


  source("~/resample/R/bootstrap2.R")
}
