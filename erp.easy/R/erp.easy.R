#' erp.easy: A user-friendly package for exploring event-related potential (ERP) data
#'
#' So you've recorded and cleaned (processed) some ERP data... Now what?
#' If you're not a programmer, or are new to ERPs, the next step may be a bit daunting, or
#' at the very least, time consuming if done by hand. The erp.easy package provides an
#' intuitive approach to exploring your processed data, without requiring a background
#' in computer programming. The erp.easy package provides three categories of functions,
#' optimized to be easy to use: loading ERP files, plotting ERP data, and analyzing ERP data.
#'
#' @section Loading functions:
#' The function \code{\link{load.data}} exists to save you the hassle of opening each individual
#' ERP file and adding a header or other identifying information to the files.  This function
#' expects data formatted with electrodes across the columns and time points as rows. Additional
#' columns (i.e., "Subject", "Stimulus", and "Time") will be added upon import to help organize your
#' data. The erp.easy package will use existing headers provided in raw data files to refer to electrodes
#' for use in functions, or will assign headers if none are present (see \code{\link{load.data}}.)
#' Single electrodes can be passed to the package functions, or several electrodes can be provided
#' (i.e., when using dense arrays) and those electrodes will be averaged together as a single electrode.
#'
#' @section Plotting functions:
#' The plotting functions \code{\link{grandaverage}}, \code{\link{individual}}, and
#' \code{\link{mosaic}} provide several ways to visualize both grand average and
#' individual data. Color-coding and labeling happens automatically for ease of use.
#'
#' @section Analysis functions:
#' The analysis functions \code{\link{m.measures}} and \code{\link{p.measures}} calculate
#' standard ERP measures such as mean amplitude, standard deviation, peak amplitude and
#' peak latency for both grand average and individual waveforms.
#'
#' @author Travis Moore
#'
#' @examples
#' library(erp.easy)
#'
#' data(ERPdata)
#'
#' grandaverage(ERPdata, electrodes = "V78")
#'
#' mosaic(ERPdata, electrodes = "V78")
#'
#' m.measures(ERPdata, electrodes = "V78", window = c(1000, 1500))
#'
#' @docType package
#' @name erp.easy
NULL
#> NULL


