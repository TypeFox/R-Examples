################################################################################
#### AUTHOR:    Arnost Komarek                                              ####
####            24/04/2014                                                  ####
####                                                                        ####
#### FILE:      confint.smoothSurvReg.R                                     ####
####                                                                        ####
#### FUNCTIONS: confint.smoothSurvReg                                       ####
####                                                                        ####
################################################################################

### ====================================================================================
### confint.smoothSurvReg: Confidence intervals for regression and scale
###                        parameters of the model fitted using the smoothSurvReg
###                        function.
### ====================================================================================

confint.smoothSurvReg <- function(object, parm, level = 0.95, method = c("pseudo-variance", "asymptotic"), ...)
{
  nParm <- nrow(object$regres)
  Parm <- rownames(object$regres)
  if (missing(parm)) parm <- Parm
  if (is.numeric(parm)){
    if (any(parm < 1) | any(parm > nParm)) stop("parm out of range.")
    parm <- Parm[parm]
  }
  if (is.character(parm)){
    if (any(!(parm %in% Parm))) stop("parm not included in the model.")
  }  

  method <- match.arg(method)
  if (level <= 0 | level >= 1) stop("level must be strictly between 0 and 1.")
  lprob <- (1 - level) / 2
  uprob <- 1 - lprob
  qq <- qnorm(uprob)

  switch(method,
    "pseudo-variance" = {
      se <- object$regres[, "Std.Error"]
    },
    "asymptotic" = {
      se <- object$regres[, "Std.Error2"]
    }
  )       
  
  low <- object$regres[, "Value"] - qq * se
  upp <- object$regres[, "Value"] + qq * se  
  names(low) <- names(upp) <- Parm
  
  if (object$estimated["common.logscale"]){
    low["Scale"] <- exp(low["Log(scale)"])
    upp["Scale"] <- exp(upp["Log(scale)"])    
  }  

  RET <- data.frame(low[parm], upp[parm])
  colnames(RET) <- paste(c(lprob, uprob) * 100, "%")

  return(RET)  
}  
