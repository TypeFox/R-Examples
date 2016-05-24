#' VWPre: Tools for Preprocessing Visual World Data.
#'
#' The VWPre package provides a set of functions for preparing Visual World data 
#' collected with SR Research Eyelink eye trackers.
#' 
#' @section Formatting functions:
#' \itemize{
#'   \item The function \code{\link{create_time_series}} returns a meaningful
#'   time columns in milliseconds.
#'   \item The function \code{\link{prep_data}} returns a data table with 
#'   correctly assigned classes for important columns.
#'   \item The function \code{\link{relabel_na}} returns a data table with
#'   samples containing 'NA' relabeled as outside any interest area.
#'   \item The function \code{\link{select_recorded_eye}} returns a data table 
#'   with data from the the recorded eye in new columns (IA_ID and IA_LABEL).
#' }
#' 
#' @section Calculation functions:
#' \itemize{
#'   \item The function \code{\link{bin_prop}} returns a downsampled data table
#'   containing proportion of looks (samples) to each interest area in a 
#'   particular window of time (bin size). 
#'   \item The function \code{\link{transform_to_elogit}} returns a data table 
#'   with proportion looks transformed to empirical logits with weights.
#'   \item The function \code{\link{create_binomial}} returns a data table with
#'   a new success/failure column for each IA which is suitable for logistic 
#'   regression. 
#' }
#' 
#' @section Fasttrack formatting function:
#' \itemize{
#'   \item The function \code{\link{fasttrack}} a meta-function that returns a 
#'   data table of processed data containing the result of the series of 
#'   necessary subroutines. This is intended for experienced users.
#' }
#' 
#' @section Utility functions:
#' \itemize{
#'   \item The function \code{\link{check_eye_recording}} returns a summary 
#'   of whether or not the dataset contains gaze data in both the Right and 
#'   Left interest area columns.  
#'   \item The function \code{\link{check_time_series}} returns the first value 
#'   in the Time column for each event. 
#'   \item The function \code{\link{check_samples_per_bin}} returns the number 
#'   of samples in each bin.
#'   \item The function \code{\link{check_samplingrate}} returns the value 
#'   corresponding to the sampling rate in the data. 
#'   \item The function \code{\link{ds_options}} returns the binning 
#'   (downsampling) options possible for the given sampling rate.
#' }
#' 
#' @section Plotting functions:
#' \itemize{
#'   \item The function \code{\link{plot_avg}} returns a plot of the grand 
#'   or conditional averages of proportion (or empirical logit) looks to each 
#'   interest area along with standard error bars.
#'   \item The function \code{\link{plot_avg_contour}} returns a contour plot of 
#'   the conditional average of proportion (or empirical logit) looks to a 
#'   given interest area over Time and a specified continuous variable. 
#' }
#' 
#' @section Diagnostic functions:
#' \itemize{
#'   \item The function \code{\link{plot_indiv_app}} opens a Shiny app for
#'   inspecting by-subject or by-item averages for all interest areas, alongside
#'   the grand average (for proportion or empirical logit looks) within a 
#'   specified time window. 
#'   \item The function \code{\link{plot_var_app}} opens a Shiny app for
#'   inspecting by-subject or by-item Z-scores with respect to the overall mean
#'   for a given interest area within a specified time window.
#' }
#' 
#' @section Notes:
#' \itemize{
#' \item The vignettes are available via \code{browseVignettes()}. 
#' \item A list of all available functions is provided in 
#' \code{help(package="VWPre")}.
#' }
#'
#' @author
#' Vincent Porretta, Jacolien van Rij, Aki-Juhani Kyröläinen, Juhani Järvikivi
#'
#' Maintainer: Vincent Porretta (\email{vincentporretta@@gmail.com})
#'
#' University of Alberta, Canada
#' 
#' @docType package
#' @name VWPre
NULL