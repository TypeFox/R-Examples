`detectAO` <-
function (object,alpha=0.05,robust=TRUE) 
{
# programmed by Kung-Sik Chan
# date: May 5, 2006
#
# This function serves to detect whether there is any AO. It implements the
# test statistic lambda_{2,t} proposed by Chang, Chen and Tiao (1988).
#
# input: 
#       object=an ARMA model
#       alpha=family significance level (5% is the default)
#             Bonferroni rule is used to control the family error rate.
#       robust=if true, the noise standard deviation is estimated by
#               mean absolute residuals times sqrt(pi/2). Otherwise,
#               it is the estimated by sqrt(sigma2) from the arima fit.
# side effects: print the test statistics of the found AO and their time indices.
#
# output:
#        a list containing the time indices, named ind, and the test statistics of the found
#        AO, named lambda2
#
#
resid=residuals(object)
piwt=ARMAtoMA(ar=-object$mod$theta,ma=-object$mod$phi,lag.max=length(resid)-1)
piwt=c(1,piwt)
omega=filter(c(0*resid[-1],rev(resid)),filter=piwt,sides=1,method='convolution')
omega=omega[!is.na(omega)]
rho2=1/cumsum(piwt^2)
omega=omega*rho2
if(robust) sigma=sqrt(pi/2)*mean(abs(resid),na.rm=TRUE) else sigma=object$sigma2^.5
lambda2T=omega/sigma/sqrt(rho2)
lambda2T=rev(lambda2T)
cutoff=qnorm(1-alpha/2/length(lambda2T))
out=abs(lambda2T)>cutoff
ind=seq(lambda2T)[out]
lambda2=lambda2T[out]
if(length(ind)!=0) print(rbind(ind, lambda2)) else print("No AO detected")
invisible(list(lambda2=lambda2, ind=ind))
}

