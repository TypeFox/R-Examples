#' @name nlbcd
#' @title Non-linearity Below Cut-off Dose
#' @usage
#' nlbcd (dosecolumn = "", targetcolumn = "",
#'        cutoffdose = 0, data = NA)
#' @description
#' This function tests non-linearity below a specified dose.
#' @details
#' The user may provide a limit below which the non-linearity of the dose-response relationship
#' is tested.  A significant result indicates that the dose-response relationship  
#' exhibits non-linearity below the user-specified cutoff dose.  NOTE: The dose-response 
#' relationship estimated by this function is not necessarily the same as that estimated by the
#' nlaad function, as the nlbcd only uses doses below the cutoff and nlaad uses all doses.
#' The user should keep this in mind in interpreting the outputs of these functions.
#' 
#' The nlaad, nlbcd, and lbcd functions are currently only intended for use on continuous outcome data.
#' @param dosecolumn   Name of dose column in input dataframe.
#' @param targetcolumn   Name of response column in input dataframe.
#' @param cutoffdose   Numeric tested cut-off dose.
#' @param data   Input dataframe.
#' @return
#' The analysis of variance table comparing the non-linear spline model with the linear model to
#' assess non-linearity across doses below the user-specified cutoff.
#' @examples
#' # Prints the F test of the difference between the spline model (output "Model 2")
#' # and the linear model (output "Model 1") as a test of nonlinearity
#' # for doses below 1.5 (i.e., all dose levels up to and including 1.49):
#' nlbcd("dose", "MF_Log", cutoffdose=1.5, data=DRdata)    
#'
#' # This produces an error, as no cutoffdose was specified:
#' \dontrun{
#' nlbcd("dose", "MF_Log", data=DRdata)
#' }
#' @export

# User can input any cutoff dose desired -- here I'm using the nogel/logel.  Can use an intermediate or truncated value to ensure only doses below cutoff get used in the following.
# This is the test of non-linearity below the NOGEL.
# Report as F with numerator df equal value in Df column and denominator df equal to Res.Df in row 2, F value and p value as indicated.
# This needs to be interpreted cautiously, as the spline is not the same as the spline from the overall analysis and can have a very different shape.

nlbcd <- function (dosecolumn = "", targetcolumn = "", cutoffdose=0, data=NA) {
    
    Below_cutoff_dose <- subset(data,data[,dosecolumn] < cutoffdose)
    targetvariablebcsubset <- Below_cutoff_dose[,targetcolumn]
    
    splinebcsubset <- mgcv::gam(targetvariablebcsubset~s(Below_cutoff_dose[,dosecolumn], k=4), data=Below_cutoff_dose)
    
    linearmodelbcsubset <- stats::lm(targetvariablebcsubset~Below_cutoff_dose[,dosecolumn], data=Below_cutoff_dose)
    
    spline_rsq <- summary(splinebcsubset)$r.sq
    linear_rsq <- summary(linearmodelbcsubset)$r.sq
    
    if ((spline_rsq - linear_rsq) > (.01*linear_rsq))
        stats::anova(linearmodelbcsubset,splinebcsubset,test="F")
    else
        return("The linear and non-linear models are not substantially different below cutoff dose.")

}
