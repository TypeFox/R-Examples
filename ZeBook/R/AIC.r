################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
################################################################################
#' @title Calculate AIC, Akaike's Information Criterion
#' @description This function calculate AIC criterion given a vector of observation, a vector of prediction and number of parameter.
#' Note that number of parameters should include variance.
#' AICcomplete is the same calculation of the AIC function of R (AICcomplete = n*log(RSS/n)+n+n*log(2*pi)+2*p, with p including variance).
#' AICshort is the calculation described in chapter 6 Working with crop model (AICshort =n*log(RSS/n)+2*p, with p including variance).
#' difference between AICcomplete and AICshort is AICcomplete-AICshort=n+n*log(2*pi)
#' As you use AIC to compare models (with different number of parameters) on a same data (with same n, number of observation), you can use AICshort or AICcomplete.
#' @param Yobs : observed values
#' @param Ypred : prediction values from the model
#' @param npar : number of parameters (should include variance that count for one supplementary parameter)
#' @return a vector with AICcomplete and AICshort
#' @export
#' @examples
#' x=c(1,2,3,4,5)
#' y=c(1.2,1.8,3.5,4.3,5.5)
#' fit = lm(y~x)
#' AIC(fit)
#' AICf(y,predict(fit),3) # 3 parameters : intercept, slope and variance
AICf<-function(Yobs,Ypred,npar){
n<-length(na.omit(Yobs))
RSS<-sum((Yobs-Ypred)^2,na.rm=TRUE)
sigma_ML=sqrt(RSS/n)
# n*log(RSS/n) + n + n*log 2pi - sum(log w)
# loglikelihood <- log((1/(2*pi*sigma_ML^2))^(n/2)) + (-RSS/(2*sigma_ML^2))
# AICcompl = -2*loglikelihood + 2*npar
AICcomplete <- n*log(RSS/n) + n + n*log(2*pi) + 2*npar
# AICde nls : -2*logLik(OLS1p) +2*2
AICshort <- n*log(RSS/n) + 2*npar
return(c(AICcomplete=AICcomplete,AICshort=AICshort))}

# end of file
