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

#' @title The state of the moving expanding window average
#'
#' @description The class holds the current state of the moving
#' expanding window (MEW) average
#'
#' @details The user should never create, update or access an instance
#' of this class themselves.  An instance of the class should be
#' created with the function \code{mewInit} and updated with the
#' functions \code{mewAccum} and \code{mewMean}.  The user can extract
#' the current value of the MEW average with the function
#' \code{mewGetMean}, and print the first six elements of the mean
#' vector to the screen with either the \code{show} or \code{print}
#' functions.
#'
#' @slot i_new (scalar integer) The index of the bin to add the
#' current sample to
#'
#' @slot i_old (scalar integer) The index of the bin to deweight
#'
#' @slot know_mean (scalar integer) flag 0: mean not known 1: mean
#' known
#'
#' @slot n_bin (scalar integer) The number of bins to use in the MEW
#' process
#'
#' @slot n_bin_use (scalar integer) The number of bins currently in
#' use
#'
#' @slot n_xx (scalar integer) The length of a vector in the sequence
#' being averaged
#'
#' @slot n_part (scalar integer) The number of samples in the bins
#' that are not being added to or deweighted
#'
#' @slot m_sample (vector integer length - n_bin) The maximum number
#' of samples allowed in each of the bins
#'
#' @slot n_sample (vector integer length - n_bin) The number of
#' samples currently in each bin
#'
#' @slot x_mean (vector double length - n_xx) The current value of the
#' MEW average (which is up-to-date only if \code{know_mean == 1})
#'
#' @slot x_sum_part (vector double length - n_xx) The sum in the bins
#' not being added to or deweighted
#'
#' @slot xx (matrix dimension - n_xx \eqn{\times} n_bin) The bin sums
#'
#' @slot ff (scalar double) The fraction of samples to retain in the
#' MEW average
#'
#' @slot ww (scalar double) The factor of increase in the number of
#' samples from one bin to the next
#'
#' @slot a_sample (scalar double) The ideal number of samples in a bin
#' (before rounding)
#'
#' @import methods
setClass(Class = "mewTyp",
         representation = representation(i_new = "integer",
             ## index of bin to add sample to (scalar)

             i_old = "integer",
             ## index of bin to deweight (if 0, don't use) (scalar)

             know_mean = "integer",
             ## flag 0: mean not known 1: mean known (scalar)

             n_bin = "integer",
             ## number of groups of samples allocated (scalar)

             n_bin_use = "integer",
             ## number of groups of samples actually defined (scalar)

             n_xx = "integer",
             ## length of each sample (scalar)

             n_part = "integer",
             ## number of groups summed (scalar)

             m_sample = "integer",
             ## max allowed samples in bin (length n_bin)

             n_sample = "integer",
             ## number of samples in each bin (length n_bin)

             x_mean = "numeric",
             ## the mean if known (length n_xx)

             x_sum_part = "numeric",
             ## sum of the part of the data which is not changing
             ## (length n_xx)

             xx = "matrix",
             ## the current bin sums (dimension n_xx rows by n_bin
             ## cols)

             ff = "numeric",
             ## fraction of samples to retain in average (typ. 1/2)
             ## (scalar)

             ww = "numeric",
             ## factor of increase in number of samples from one bin
             ## to the next (scalar)

             a_sample = "numeric"
             ## ideal number of samples in a bin (before rounding)
             ## (scalar)
             ))
