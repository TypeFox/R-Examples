estSSR <-
function(x, y, family="normal", twoside=TRUE, type="RG", alpha=0.05, B=2000)
{
# x: sample from r.v. X
# y: sample from r.v. Y
# family: distribution of X and Y
# twoside: TRUE for two-sided CI; FALSE for one-sided CI
# type: type of CI to be built
# alpha: 1-alpha is the nominal confidence level of the CI to be built
# B: number of bootstrap replicates (for the bootstrap CI only)
if(family=="normal")
{
#sample quantities
mux<-mean(x)
muy<-mean(y)
n<-length(x)
m<-length(y)
S2x<-var(x)
S2y<-var(y)
S2<-S2x+S2y
# MLE
hatd<-(mean(x)-mean(y))/sqrt(S2x*(n-1)/n+S2y*(m-1)/m)
hatR<-pnorm(hatd)
# DOWNTON
hatd2<-(mean(x)-mean(y))/sqrt(S2x+S2y)
hatR2<-pnorm(hatd2)
# construction of confidence intervals
# Reiser-Guttman
if(type=="RG")
{
delta<-(mux-muy)/sqrt(S2)
M<-(S2x+S2y)/(S2x/n+S2y/m)
f<-(S2x+S2y)^2/(S2x^2/(n-1)+S2y^2/(m-1))
if(twoside)
{
delta1<-delta-(1/M+delta^2/(2*f))^0.5*qnorm(1-alpha/2)
delta2<-delta+(1/M+delta^2/(2*f))^0.5*qnorm(1-alpha/2)
R1<-pnorm(delta1)
R2<-pnorm(delta2)
CI<-c(R1,R2)
}
else
{
delta2<-delta-(1/M+delta^2/(2*f))^0.5*qnorm(1-alpha)
CI<-pnorm(delta2)
}
}
# Guo-Krishnamoorthy
else if(type=="GK")
{
delta<-(mux-muy)/sqrt(S2)
hatq1<-S2x*(m-3)/(S2y*(m-1))
hatq2<-S2y*(n-3)/(S2x*(n-1))
hatm1<-n*(1+hatq1)/(hatq1+n/m)
hatm2<-m*(1+hatq2)/(hatq2+m/n)
hatf1<-(n-1)*(hatq1+1)^2/(hatq1^2+(n-1)/(m-1))
hatf2<-(m-1)*(hatq2+1)^2/(hatq2^2+(m-1)/(n-1))
delta1l<-gkf(1-alpha,sqrt(hatm1)*delta,hatf1)/sqrt(hatm1)
delta2l<-gkf(1-alpha,sqrt(hatm2)*delta,hatf2)/sqrt(hatm2)
deltagk<-min(delta1l,delta2l)
CI<-pnorm(deltagk)
}
# large-sample asymptotic normal
else if (type=="AN")
{
hatVhatd<-1/(S2)*(S2x/n+S2y/m+0.5/S2^2*(S2x^2/n+S2y^2/m)*(mean(x)-mean(y))^2)
hatVhatd*1/(2*pi*S2)*exp(-(mean(x)-mean(y))^2/(S2))
if(twoside)
{
qdt<-c(hatd-qnorm(1-alpha/2)*sqrt(hatVhatd),hatd+qnorm(1-alpha/2)*sqrt(hatVhatd))
CI<-pnorm(qdt)
}
else
{
qdt<-hatd-qnorm(1-alpha)*sqrt(hatVhatd)
CI<-pnorm(qdt)
}
}
# parametric percentile bootstrap
else if (type=="B")
{
Rb<-numeric(B)
for(j in 1:B)
{
xb<-rnorm(n,mux,sqrt(S2x*(n-1)/n))
yb<-rnorm(m,muy,sqrt(S2y*(m-1)/m))
muxb<-mean(xb)
muyb<-mean(yb)
sigmaxb<-sqrt(var(xb)*(n-1)/n)
sigmayb<-sqrt(var(yb)*(m-1)/m)
Rb[j]=pnorm((muxb-muyb)/sqrt(sigmaxb^2+sigmayb^2))
}
if(twoside)
{
q<-quantile(Rb,c(alpha/2,1-alpha/2))
CI<-q
}
else
{
q<-quantile(Rb,alpha)
CI<-q
}
}
# variance stabilizing transformation - logit
else if (type=="LOGIT")
{
hatVhatRe<-1/(2*pi*S2)*exp(-(mean(x)-mean(y))^2/(S2))*(S2x/n+S2y/m+0.5/S2*(S2x^2/S2/n+S2y^2/S2/m)*(mean(x)-mean(y))^2)
logitR<-log(hatR/(1-hatR))
if(twoside)
{
qtlogit<-c(logitR-qnorm(1-alpha/2)*sqrt(hatVhatRe)/(hatR*(1-hatR)),logitR+qnorm(1-alpha/2)*sqrt(hatVhatRe)/(hatR*(1-hatR)))
CI<-exp(qtlogit)/(1+exp(qtlogit))
}
else
{
qtlogit<-logitR-qnorm(1-alpha)*sqrt(hatVhatRe)/(hatR*(1-hatR))
CI<-exp(qtlogit)/(1+exp(qtlogit))
}
}
# variance stabilizing transformation - arcsin
else if (type=="ARCSIN")
{
hatVhatRe<-1/(2*pi*S2)*exp(-(mean(x)-mean(y))^2/(S2))*(S2x/n+S2y/m+0.5/S2*(S2x^2/S2/n+S2y^2/S2/m)*(mean(x)-mean(y))^2)
logitsin<-asin(sqrt(hatR))
if(twoside)
{
qtsin<-c(max(0,logitsin-qnorm(1-alpha/2)*sqrt(hatVhatRe/(4*hatR*(1-hatR)))),min(pi/2,logitsin+qnorm(1-alpha/2)*sqrt(hatVhatRe/(4*hatR*(1-hatR)))))
CI<-(sin(qtsin))^2
}
else
{
qtsin<-logitsin-qnorm(1-alpha)*sqrt(hatVhatRe/(4*hatR*(1-hatR)))
CI<-(sin(qtsin))^2
}
}
}
# return results as a list
return(list(ML_est=hatR,Downton_est=hatR2,CI=CI,confidence_level=1-alpha))
}

