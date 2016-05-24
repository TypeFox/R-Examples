## This is a Translation of Zachary Levine's mewAvg.f90 to R All of
## Zachary's comments are included in this code.  With the exception
## of the roxygen2 documentation, unless a comment is specifically
## noted to be Adam's, it is Zachary's

## Zachary Levine 23 May 2013 - 4 June 2013
## The purpose of this module is to implement a particular averaging
## scheme which allows convergence in stochastic optimization.  The
## background is described in JC Spall, "Intro. to Stochastic Search
## and Optimization" Wiley, 2003, Chap. 4.

## The idea is to obtain the average from the following sum (for N
## even):

## \bar X = \lim_{N->\infty} {2/N} \sum_{i=(N/2)+1}^N X_i

## That is, we have a moving, expanding average, where the first half
## of the samples are discarded.  The first half are discarded because
## they are obtained under conditions which are different than those
## of the converged parameters. (The "half" is parameterized in the
## implementation.)

## In order to use a fixed amount of storage as N->\infty, we will
## have a fixed number of bins (nBin) which are partial sums of the
## series.  The number of samples in each bin increases exponentially
## (by a factor of ww, rounded to an integer).  The oldest bin is
## phased out as the newest bin is filled.

## To avoid keeping track of many shapes, the X_i is taken to be a 1D
## array.

## At the begining, only one sample is stored per bin until all bins
## have at least one sample.  At the very beginning, the mean is set
## to 0.

## Usage
##   loop over independent uses
##      call mewInit
##      loop over sample acquisition and use of mean
##         call mewAccum (when new data exists)
##         call mewMean  (whenever desired)
##      call mewFinal (optional - space reuseable in any case)

## Adam Comment
## This function initializes an instance of the class mewTyp
## An object of class mewTyp should never be initialized by the user
## with the "new" function
## If mewInit is successful, a correctly initialized instance of the
## class mewTyp will be returned; otherwise, an error is raised

#' @title Create an S4 object of class \code{mewTyp}
#'
#' @description Call this function to create an S4 object of class
#' \code{mewTyp}.
#'
#' @details If it is necessary to directly call mewAccum and mewMean
#' an S4 object of class \code{mewTyp} should be created first using
#' this function.  The user should never create an S4 object of class
#' \code{mewTyp} using the \code{new} function provided by the
#' \code{methods} package.
#'
#' @param n_bin (scalar integer) The fixed number of bins to use to
#' define the moving expanding window
#'
#' @param n_xx (scalar integer) The length of each vector in the
#' sequence to be averaged
#'
#' @param ff (scalar double) The fraction of the samples to included
#' in each window
#'
#' @return An initialized instance of the class \code{mewTyp}
#'
#' @examples
#' av <- mewInit(n_bin = 4, n_xx = 2, ff = 0.5)
#'
#' @export
mewInit <- function (n_bin, n_xx, ff) {

  av <- new("mewTyp")

  av@n_bin <- as.integer(n_bin)
  av@n_xx <- as.integer(n_xx)
  av@ff <- as.double(ff)

  ## checking input for errors
  if (av@n_bin < 3) {

    stop(paste0("mewInit: need 3<=n_bin but was ", av@n_bin))
  }

  if ((av@ff <= 0.0) || (av@ff >= 1.0)) {

    stop(paste0("mewInit: need 0<ff<1, but was ", av@ff))
  }

  if (av@n_xx <= 0) {

    stop(paste0("mewInit: need n_xx>0 but was ", av@n_xx))
  }

  ## increment factor
  av@ww <- (1.0 - av@ff)^(-1.0/av@n_bin)

  ## ideally this many samples
  av@a_sample <- 1.0

  ## no bins in use so far
  av@n_bin_use <- as.integer(0)

  ## first bin to write minus 1
  av@i_new <- as.integer(0)

  ## first bin to overwrite
  av@i_old <- as.integer(1)

  ## mean known
  av@know_mean <- as.integer(1)

  ## Adam comment
  ## allocate memory for x_mean (of length n_xx)
  ## initialize x_mean
  av@x_mean <- double(av@n_xx)

  ## Adam comment
  ## allocate memory for x_sum_part (of length n_xx)
  ## initialize x_sum_part
  av@x_sum_part <- double(av@n_xx)

  ## Adam comment
  ## allocate memory for xx (dimension n_xx rows by n_bins cols)
  ## initialize xx
  av@xx <- matrix(data = 0,
                  nrow = av@n_xx,
                  ncol = av@n_bin)

  ## Adam comment
  ## allocate memory for m_sample (length n_bin)
  ## initialize m_sample
  av@m_sample <- as.integer(rep(x = 1, times = av@n_bin))

  ## Adam comment
  ## allocate memory for n_sample (length n_bin)
  ## initialize n_sample
  av@n_sample <- integer(av@n_bin)

  av@n_part <- as.integer(0)

  return(av)
}
