#' @name nlaad
#' @title Non-linearity Across All Doses
#' @usage
#' nlaad(dosecolumn = "", targetcolumn = "", data = NA)
#' @description
#' This function determines whether a non-linear spline model of the dose-response
#' relationship differs significantly from a linear model across all doses.
#' @details
#' The non-linear spline model (output "Model 2") is compared to the linear model (output "Model 1") 
#' using an anova F-test.  If the spline model fits the data signficantly better,
#' the F will be large and the associated p value will be significant.
#' 
#' The nlaad, nlbcd, and lbcd functions are currently only intended for use on continuous outcome data.
#' @param dosecolumn   Name of dose column in dataframe.
#' @param targetcolumn   Name of response column in dataframe.
#' @param data   Input dataframe.
#' @return
#' The analysis of variance table comparing the non-linear spline model with the linear model to
#' assess non-linearity across all doses.
#' @examples
#' # Prints the F test of the difference between the spline model (output "Model 2")
#' # and the linear model (output "Model 1") as a test of nonlinearity
#' # across all doses:
#' nlaad("dose", "MF_Log", data=DRdata)
#' @export

# Report as F with numerator df equal value in Df column and denominator df equal to Res.Df in row 2, F value and p value as indicated.

nlaad <- function (dosecolumn="", targetcolumn="", data=NA) {
    
    splinemodel <- mgcv::gam(data[,targetcolumn]~s(data[,dosecolumn], k=4), data=data)
    linearmodel <- stats::lm(data[,targetcolumn]~data[,dosecolumn], data=data)
    
    spline_rsq <- summary(splinemodel)$r.sq
    linear_rsq <- summary(linearmodel)$r.sq
    
    if ((spline_rsq - linear_rsq) > (.01*linear_rsq))
        stats::anova(linearmodel,splinemodel,test="F")
    else
        return("The linear and non-linear models are not substantially different.")
}


