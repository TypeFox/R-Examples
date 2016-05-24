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

## ###################################################################
## Writes the mean using a windowed moving average.
## The oldest data is reduced in weight linearly as the new data comes
## in.

#' @title Update the moving expanding window average
#'
#' @description When desired, the \code{x_mean} slot in an S4 object
#' of class \code{mewTyp} may be updated to contain the correct moving
#' expanding window (MEW) average (it is not updated by the function
#' \code{mewAccum} to save computation).  If the slot \code{know_mean}
#' is unity, the slot \code{x_mean} is up-to-date; otherwise; it is
#' not.
#'
#' @param av (class mewTyp) the current state of the MEW average
#'
#' @return the updated instance of the argument \code{av}
#'
#' @examples
#' n_iter <- 100
#'
#' i_to_print <- 10
#'
#' results <- matrix(data = double(2*n_iter/i_to_print),
#'                   nrow = n_iter/i_to_print,
#'                   ncol = 2)
#'
#' av <- mewInit(n_bin = 4, n_xx = 2, ff = 0.5)
#'
#' for (i in 1:n_iter) {
#'
#'   value <- runif(n=2)
#'   value[1] <- ((cos(value[1]*2*pi))^2)*(1 - exp(-0.01*i))
#'   value[2] <- (-((sin(value[2]*2*pi))^2))*(1 - exp(-0.01*i))
#'
#'   av <- mewAccum(xx = value, av = av)
#'
#'   if (i%%i_to_print == 0) {
#'
#'     av <- mewMean(av)
#'     show(av)
#'     results[i/i_to_print, ] <- mewGetMean(av)
#'   }
#' }
#'
#' ## plot the results
#'
#' plot(c(1, (n_iter/i_to_print)),
#'      c(min(results), max(results)),
#'      type = "n")
#' points(1:(n_iter/i_to_print), results[, 1])
#' points(1:(n_iter/i_to_print), results[, 2])
#'
#' ## Now, a larger example, and we pause part way through to assess
#' ## convergence
#'
#' n_iter <- 1000
#' av <- mewInit(n_bin = 4, n_xx = 5000, ff = 0.5)
#' for (i in 1:n_iter) {
#'
#'   new_samp <- runif(n = 5000)
#'   av <- mewAccum(xx = new_samp, av = av)
#' }
#'
#' av <- mewMean(av = av)
#'
#' ## of course each element of the mean sould converge to 0.5.  After
#' ## 1000 iterations, the first six elements of the mean vector are
#' show(av)
#'
#' ## run another 1000 iterations
#' for (i in 1:1000) {
#'
#'   new_samp <- runif(n = 5000)
#'   av <- mewAccum(xx = new_samp, av = av)
#' }
#'
#' av <- mewMean(av)
#'
#' ## check the mean of the first six elements again
#' show(av)
#'
#' @useDynLib mewAvg
#'
#' @export
mewMean <- function(av) {

  ## Adam comment
  ## checking the argument for class mewTyp
  if (class(av)[1] != "mewTyp") {

    stop("mewMean: the argument should be of class mewTyp")
  }

  if (av@know_mean == as.integer(1)) {
    ## if the mean is known, no action is needed

    return(av)
  }

  i_new <- av@i_new
  i_old <- av@i_old

  av@know_mean <- as.integer(1)

  fract <- (av@m_sample[i_new] + 1.0 - av@n_sample[i_new])/
    (av@m_sample[i_new] + 1.0)

  .Call("meanCalc",
        av@x_mean,
        av@xx,
        av@x_sum_part,
        fract,
        av@n_xx,
        as.integer(i_old - 1), ## C indexes start at zero (Adam)
        as.integer(i_new - 1), ## C indexes start at zero (Adam)
        av@n_sample,
        av@n_part)

  return(av)
}
