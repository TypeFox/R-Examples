
VSPKnown=function(N, n, k, Pa=seq(0, 1, 0.001), Plots=TRUE)
{
zpa = qnorm(Pa)
zp = k+(zpa/sqrt(n))
p=1-pnorm(zp)
OC = Pa
AOQ=(N-n)*p*OC/N
ATI= n*OC+N*(1-OC)
results = list(p=p, OC=OC, n=rep(n,length(p)),  AOQ=AOQ, ATI=ATI)
class(results)="AccSampPlan"
if(Plots){
par(mfrow=c(2,2))
plot(results)
}
return(results)
}

VSPUnknown=function(N, n, k, Pa=seq(0, 1, .001), Plots=TRUE)
{
zpa = qnorm(Pa)
k1 = sqrt(1+(k*k)/2)
zp = k+(k1*zpa/sqrt(n))
p=1-pnorm(zp)
OC = Pa
AOQ=(N-n)*p*OC/N
ATI= n*OC+N*(1-OC)
results = list(p=p, OC=OC, n=rep(n,length(p)),  AOQ=AOQ, ATI=ATI)
class(results)="AccSampPlan"
if(Plots){
par(mfrow=c(2,2))
plot(results)
}
return(results)
}




VSPDesign =function(AQL, alpha, LQL, beta)
{
zp1=qnorm(1-AQL)
zp2=qnorm(1-LQL)
zpa1=qnorm(1-alpha)
zpa2=qnorm(1-beta)
k =(zp2*zpa1+zp1*zpa2)/(zpa1+zpa2)
n =(zpa1+zpa2)/(zp1-zp2)
n = n*n
n = round(n)
n.unknown=n*(1+(k*k/2))
return(data.frame(k, n, n.unknown))
}

