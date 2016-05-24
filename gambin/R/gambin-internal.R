.est_confint <-
function(est_alpha, est_likelihood, mydata, level)
{
  conf_logLik <- function(alpha, mydata, est_likelihood, level)
  {
    return(abs(.logLik_gamBin(alpha, mydata) - est_likelihood - exp(-qchisq(level,1)/2))) # the likelihood ratio test for confidence intervals (from "Beyond Traditional Statistical Measures")
  }
  lower <- optimise(conf_logLik, interval = c(0,est_alpha), mydata = mydata, est_likelihood = est_likelihood, level = level)$minimum
  higher <- optimise(conf_logLik, interval = c(est_alpha, 30), mydata = mydata, est_likelihood = est_likelihood, level = level)$minimum
  return(c(lower, higher))
}
.logLik_gamBin <-
function(alpha, mydata, maxoct) #the maxoct argument can be removed if we do not estimate maxoct
{
  if(missing(maxoct)) maxoct <- max(mydata$octave)
  dgamb <- dgambin(alpha, maxoct)
  exponent <- mydata$species # this line and the next can be removed if we do not estimate maxoctave
  if(length(exponent) < length(dgamb)) exponent[(length(exponent)+1):length(dgamb)] <- 0
  lik_dist <- dgamb^exponent
  logLik <- sum(-log(lik_dist))
  logLik
}
.sample_abundances <-
function(abundances, individuals)
{
  if(!is.numeric(abundances)) stop("abundances must be numeric")
  species <- as.factor(1:length(abundances))
  samplevector <- rep(species, abundances)
  samp <- sample(samplevector, individuals, replace = F)
  return(table(samp))
}
