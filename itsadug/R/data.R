#' Simulated time series data.
#'
#' A dataset containing the sine wave data with random noise added. 
#'
#' @format A data frame with 75600 rows and 6 variables:
#' \describe{
#'   \item{\code{Group}}{Age group of participants: Adults or Children.}
#'   \item{\code{Time}}{Time, time measure from start of each time series.}
#'   \item{\code{Trial}}{Trial in the experiment, centered around zero.}
#'   \item{\code{Condition}}{Continuous variable, ranging from -1 to 4. 
#' For example, stimulus onset asynchrony.}
#'   \item{\code{Subject}}{Code for individual participants.}
#'   \item{\code{Y}}{Time series measure. 
#' Similar to pupil size, sensor position, or voltage.}
#' }
#' @author Jacolien van Rij
"simdat"
#' Raw EEG data, single trial, 50Hz.
#'
#' A dataset containing a single EEG trial.
#'
#' @format A data frame with 1504 rows and 5 variables:
#' \describe{
#'   \item{\code{Electrode}}{Electrode that recorded the EEG.}
#'   \item{\code{Time}}{Time, time measure from onset of the stimulus.}
#'   \item{\code{Ampl}}{EEG amplitude, recorded by 32 electrodes.}
#'   \item{\code{X}}{Approximation of electrode position, relative to Cz. 
#' Left is negative.}
#'   \item{\code{Y}}{Approximation of electrode position, relative to Cz. 
#' Back is negative.}
#' }
#' @author Jacolien van Rij
"eeg"





