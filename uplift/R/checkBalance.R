######################################################################
# Wrapper function for xBalance{RItools}
# Test conditional independence of the treatment variable and the 
# covariates.
######################################################################

checkBalance <- function(formula, data, report = "all", ...) {
  
  xBalance(fmla = formula, data = data, report = report, ...)
  
}