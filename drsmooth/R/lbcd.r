#' @name lbcd
#' @title Linearity Below Cut-off Dose
#' @usage
#' lbcd(dosecolumn = "", targetcolumn = "", cutoffdose = 0, data = NA)
#' @description
#' This function tests linearity below a specified dose.
#' @details
#' The user may provide a limit below which the linearity of the dose-response relationship
#' is tested.  A significant result indicates that the slope is non-zero below 
#' the user-specified cutoff dose.
#' 
#' The nlaad, nlbcd, and lbcd functions are currently only intended for use on continuous outcome data.
#' @param dosecolumn   Name of dose column in input dataframe.
#' @param targetcolumn   Name of response column in input dataframe.
#' @param cutoffdose   Cut-off dose (numeric).
#' @param data   Input dataframe.
#' @return
#' A summary table showing the intercept and slope coefficients for the linear function below
#' the user-specified dose, along with standard errors and significance tests.
#' @examples
#' # Conducts a linear regression for all doses below 1.5
#' # (i.e., all dose levels up to and including 1.4929).
#' # The significance test on the dose coefficient is the test of non-zero linear slope:
#' lbcd("dose", "MF_Log", cutoffdose=1.5, data=DRdata)
#'
#' # This produces an error, as no cutoffdose was specified:
#' \dontrun{
#' lbcd("dose", "MF_Log", data=DRdata)
#' }
#' @export

lbcd <- function (dosecolumn = "", targetcolumn = "", cutoffdose=0, data=NA) {
    
    Dose_below_cutoff <- subset(data,data[,dosecolumn] < cutoffdose)
    targetvariablebcsubset <- Dose_below_cutoff[,targetcolumn]
    
    linearmodelbcsubset <- stats::lm(targetvariablebcsubset~Dose_below_cutoff[,dosecolumn], data=Dose_below_cutoff)

    summary(linearmodelbcsubset)
}
