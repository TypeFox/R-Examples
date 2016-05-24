#' @title Dose-response Modeling with Smoothing Splines
#' @usage
#' smooth(dosecolumn = "", 
#'        targetcolumn = "", 
#'        data = NA)
#' @description
#' Replaced by drsmooth::drsmooth.
#' @details
#' drsmooth::smooth has been replaced by drsmooth::drsmooth but a shell retained
#' in v.1.9.0 to prompt the switch.  v.2.0.0 will remove drsmooth::smooth altogether.
#' Please see help for drsmooth, where additional detail on the new functionality
#' is also available.
#' @param dosecolumn   Name of dose column of interest in dataframe.
#' @param targetcolumn  Name of response column of interest in dataframe.
#' @param data   Input dataframe.
#' @export

smooth <- function (dosecolumn="", 
                      targetcolumn="", 
                      data=NA) {
  .Deprecated("drsmooth")
  warning("Please pardon this non-backward compatible change, to prepare for a correction for the conflict with stats::smooth.", 
          call. = FALSE)
}