#############################################
#### AUTHOR:    Arnost Komarek           ####
####            01/05/2004               ####
####                                     ####
#### FILE:      summary.smoothSurvReg.R  ####
####                                     ####
#### FUNCTIONS: summary.smoothSurvReg    ####
#############################################

### ===========================================================================
### summary.smoothSurvReg: Print summary for objects of class 'smoothSurvReg'
### ===========================================================================
## object ..... object of class 'smoothSurvReg'
## spline ..... T/F, do I want to print an information concerning the fitted spline?
## digits ..... # of printed digits
## ... ........ other arguments passed to 'print' function
summary.smoothSurvReg <- function(object, spline, digits = min(options()$digits, 4), ...)
{
   print.smoothSurvReg(object, spline, digits, ...)
}

