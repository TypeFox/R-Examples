# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

samp.bootstrap <- function(n, R, size = n - reduceSize, reduceSize = 0)
{
  # Generate R random samples from 1:n
  #
  # Args:
  #   n: integer, draw values with replacement from 1:n.
  #   R: number of permutations.
  #   size: if specified, then samples of size 'size'.
  #   reduceSize: an alternate way of specifying size.
  #
  # Return:
  #   a matrix with R columns and size rows, of samples drawn with replacement.
  matrix(ceiling(n * runif(size * R)),
         nrow = size, ncol = R)
}

samp.bootknife <- function(n, R, size = n - reduceSize, reduceSize = 0)
{
  # Generate R random samples from 1:n, each from a jackknife sample.
  #
  # Args:
  #   n: integer, draw values with replacement from 1:n.
  #   R: number of permutations.
  #   size: if specified, then samples of size 'size'.
  #   reduceSize: an alternate way of specifying size.
  #
  # Return:
  #   a matrix with R columns and size rows, of samples drawn with replacement.
  omit <- ceiling(n * runif(R)) # randomly omit observations
  result <- matrix(ceiling((n-1) * runif(size * R)),
                   nrow = size, ncol = R)
  # Where an observation equals a value that should be omitted, change to n.
  result[result == rep(omit, each = size)] <- n
  result
}

samp.permute <- function(n, R, size = n - reduceSize, reduceSize = 0,
                         groupSizes = NULL, returnGroup = NULL)
{
  # Generate R random permutations from 1:n
  #
  # Args:
  #   n: integer, draw values without replacement from 1:n.
  #   R: number of permutations.
  #   size: if specified, then partial permutations, of size values from 1:n.
  #         size is ignored if groupSizes and returnGroup are given..
  #   reduceSize: an alternate way of specifying size.
  #   groupSizes: vector of positive integers, that sums to n.
  #               Let k = its length.
  #               For a k-sample problem, specify groupSizes and returnGroup.
  #   returnGroup: integer from 1:k, indicating which group to return.
  #
  # Return:
  #   a matrix with R columns and size or groupSizes[returnGroup] rows,
  #   where each column is a (partial) permutation of 1:n.
  #
  # groupSizes and returnGroup are useful in k-sample problems. By
  # resetting .Random.seed to the same value before each of k calls,
  # and letting returnGroup be values from 1:k, the selections for
  # different groups can be obtained in successive calls.

  if(xor(is.null(groupSizes), is.null(returnGroup)))
    stop("groupSizes and returnGroup must both be specified, or neither.")
  if(!is.null(returnGroup)) {
    # Generate full permutations, but only return part
    cumSizes <- c(0, cumsum(groupSizes))
    indices <- seq(from = cumSizes[returnGroup] +1,
                   to =   cumSizes[returnGroup+1], by = 1)
    result <- matrix(integer(R * groupSizes[returnGroup]), ncol = R)
    for(j in 1:R)
      result[, j] <- sample(n)[indices]
    return(result)
  }
  stopifnot(size <= n, size > 0)
  result <- matrix(integer(R * size), ncol = R)
  for(j in 1:R)
    result[, j] <-  sample(n)[1:size]
  # That is not efficient if size < n, but does return the same values
  # as when size = n.
  return(result)
}
# TODO: more efficient sampling, handle the loop inside C.



if(FALSE) {
  source("~/resample/R/samplers.R")

  set.seed(0)
  i1 <- samp.bootstrap(10, 12)
  set.seed(0)
  i2 <- matrix(ceiling(10 * runif(120)), 10)
  all.equal(i1, i2)

  set.seed(0)
  i3 <- matrix(sample(10, size = 120, replace = TRUE), 10)
  all.equal(i1, i3)

  set.seed(0)
  i4 <- samp.bootknife(10, 12)

  var(colMeans(samp.bootstrap(10, 10^4)))
  var(colMeans(samp.bootstrap(10, 10^4, reduceSize = 1)))
  var(colMeans(samp.bootknife(10, 10^4)))

  set.seed(0); j1 <- samp.permute(10, 12)
  set.seed(0); j2 <- samp.permute(10, 12, size = 6)
  set.seed(0); j3 <- samp.permute(10, 12, reduceSize = 4)
  set.seed(0); j4 <- samp.permute(10, 12, groupSizes = c(6, 4), returnGroup = 1)
  set.seed(0); j5 <- samp.permute(10, 12, groupSizes = c(6, 4), returnGroup = 2)
  all.equal(j1[1:6, ], j2)
  all.equal(j2, j3)
  all.equal(j1[1:6, ], j4)
  all.equal(j1[7:10, ], j5)
}
