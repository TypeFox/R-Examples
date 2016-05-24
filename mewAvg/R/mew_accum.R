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

## ##################################################################/
##    Accumulate samples for average
## (1) not all bins are in use yet; assign the new sample to the next
##     bin
## (2) all bins are in use
##     (2a) accumulate in existing bin
## (2b) eliminate an old bin, and start a new one.  Also, the partial
##      sum of all the not-old-and-not-new bins is found.
##     (2c) error
## (3) error

#' @title Update the class \code{mewTyp}
#'
#' @description Update an S4 object of class \code{mewTyp} with a new
#' data point
#'
#' @details If \code{av} is an S4 object of class \code{mewTyp} that
#' contains the current state of the MEW average and \code{xx} is a
#' new vector of data, the function \code{mewAccum} updates the MEW
#' average with \code{xx}.
#'
#' @param xx (vector double) The vector of data with which to update
#' the MEW aveage
#'
#' @param av (class mewTyp) The current state of the MEW average
#'
#' @return The updated instance of \code{av}
#'
#' @examples
#' n_iter <- 1000
#'
#' av <- mewInit(n_bin = 4, n_xx = 1, ff = 0.5)
#'
#' for (i in 1:n_iter) {
#'
#'   value <- runif(n=2)
#'   value[1] <- ((cos(value[1]*2*pi))^2)*(1 - exp(-0.01*i))
#'   value[2] <- (-((sin(value[2]*2*pi))^2))*(1 - exp(-0.01*i))
#'   value <- as.double(value)
#'
#'   av <- mewAccum(xx = value, av = av)
#' }
#'
#' @export
#'
#' @useDynLib mewAvg
mewAccum <- function (xx, av) {

  ## Adam comment
  ## checking the first argument for type double
  if (!is.double(xx)) {

    stop("mewAccum: the first argument should be of type double")
  }

  ## Adam comment
  ## checking the second argument for class mewTyp
  if (class(av)[1] != "mewTyp") {

    stop("mewAccum: the second argument should be of class mewTyp")
  }

  ## Adam comment
  ## unpack to simplify syntax
  i_new <- av@i_new
  i_old <- av@i_old

  i_not_new_not_old <- integer(av@n_bin - 2)

  if (av@n_bin_use < av@n_bin) {
    ## (1)

    i_new <- as.integer(i_new + 1)
    av@n_bin_use <- as.integer(av@n_bin_use + 1)

    .Call("replaceCol",
          av@xx,
          xx,
          as.integer(i_new - 1), ## C index starts at zero
          av@n_xx)

    av@n_sample[i_new] <- as.integer(1)

    if (av@n_bin_use == as.integer(1)) {

      .Call("assignLongVec",
            av@x_mean,
            xx,
            av@n_xx)
    } else {

      .Call("binNotFullMean",
            av@x_mean,
            xx,
            av@n_xx,
            av@n_bin_use)
    }
  } else if (av@n_bin_use == av@n_bin) {
    ## (2)

    av@know_mean <- as.integer(0)

    if (av@n_sample[i_new] < av@m_sample[i_new]) {
      ## (2a)

      av@n_sample[i_new] <- as.integer(av@n_sample[i_new] + 1)

      .Call("addToBin",
            av@xx,
            xx,
            as.integer(i_new - 1), ## C index starts at zero
            av@n_xx)
    } else if (av@n_sample[i_new] == av@m_sample[i_new]) {
      ## (2b)

      i_new <- i_old

      if (i_old < (av@n_bin)) {

        i_old <- as.integer(i_old + 1)
      } else {

        i_old <- as.integer(1)
      }

      .Call("replaceCol",
            av@xx,
            xx,
            as.integer(i_new - 1), ## C index starts at zero
            av@n_xx)

      av@n_sample[i_new] <- as.integer(1)
      av@a_sample <- av@a_sample*av@ww
      av@m_sample[i_new] <- as.integer(round(av@a_sample))

      j <- as.integer(0)
      for (i in 1:av@n_bin) {

        if ((i == i_old) || (i == i_new)) {

          next
        }

        j <- as.integer(j + 1)
        i_not_new_not_old[j] <- as.integer(i)
      }

      av@n_part <- as.integer(sum(av@n_sample[i_not_new_not_old]))

      .Call("addxSumPart",
            av@x_sum_part,
            av@xx,
            i_not_new_not_old,
            av@n_xx,
            av@n_bin)
    } else {
      ## (2c)

      cat(paste0("mewAccum: likely programming error i_new is ",
                 i_new, "\n"))
      cat(paste0("mewAccum: n_sample[i_new] is ",
                 av@n_sample[i_new], "\n"))
      cat("mewAccum: need <= \n")
      cat(paste0("mewAccum: m_sample[i_new] which is ",
                 av@m_sample[i_new], "\n\n"))
      stop()
    }
  } else {
    ## (3)

    cat("mewAccum: likely programming error\n")
    cat(paste0("mewAccum: n_bin_use is ",
               av@n_bin_use, "\n"))
    cat("mewAccum: need <= \n")
    cat(paste0("mewAccum: n_bin which is ",
               av@n_bin, "\n\n"))
    stop()
  }

  av@i_new <- as.integer(i_new)
  av@i_old <- as.integer(i_old)

  return(av)
}

