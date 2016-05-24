# this function provides CI for a GAMLSS model
# it relies on the vcov() working which I hope is the case for parameteric models
# TO DO 
# i) worth inversigating the function profile.glm() in MASS MASS profile.glm
# since it is the function  used for glm models and used in confint.glm()
# ------------------------------------------------------------------------------
confint.gamlss  <- function (object, parm, level = 0.95, 
                         what = c("mu", "sigma", "nu", "tau"), 
                         parameter = NULL, robust = FALSE, ...) 
{
#local function from stats	
format.perc <- function (probs, digits) 
paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), 
    "%")	
#--------------------------	
 what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
  cf <- coef(object, what)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             pct))
  ses <- sqrt(diag(vcov(object, robust=robust)))[parm] # this needs checking
  ci[] <- cf[parm] + ses %o% fac
  ci
}