

SeqDesignBinomial=function(N=NULL, AQL, alpha, LQL, beta, Plots=TRUE){
a = log((1-beta)/alpha)
b = log((1-alpha)/beta)
g1= log(LQL/AQL)
g2 = log((1-AQL)/(1-LQL))
G=g1+g2
h1= b/G
h2 = a/G
s = g2/G

h = seq(- 4*h1, 5*h2, 0.01)
p =(1-((1-LQL)/(1-AQL))^h )/(((LQL/AQL)^h)-(((1-LQL)/(1-AQL))^h))
L = round(2*h1/s)
k = seq(1, L, 1)
accept = s*k-h1
reject = s*k+h2

results=list(N=N, AQL=AQL, alpha=alpha, LQL=LQL, beta=beta, h1=h1, h2=h2, s=s, accept=accept, reject=reject, k=k,p=p, h=h)
class(results)="SeqSampPlan"
if(Plots){
plot(results)
}
return(results)
}


SequentialBinomial=function(x, Plots=TRUE){
k1=((1-x$beta)/x$alpha)^x$h
k2 =(x$beta/(1-x$alpha))^x$h
OC =(k1-1)/(k1-k2)
AOQ= x$p*OC

k5 = OC* log(x$beta/(1-x$alpha)) 
k3=k5+(1-OC)* log((1-x$beta)/x$alpha)
k4 = x$p*log(x$LQL/x$AQL)+(1-x$p)*log((1-x$LQL)/(1-x$AQL))
ASN = k3/k4
ATI=k5/k4 + (1-OC)*x$N
results = list(p=x$p, OC=OC, ASN=ASN, AOQ=AOQ, ATI=ATI)
class(results)="AccSampPlan"
if(Plots){
par(mfrow=c(2,2))
plot(x$p, OC, type="l", 
ylab="Probability of Acceptance", 
xlab="Fraction Nonconforming p")
title("OC Curve")

plot(x$p, ASN, type="l", 
ylab="Average sample size", 
xlab="Fraction Nonconforming p")
title(paste("maximum ASN = ", formatC(max(ASN))))

plot(x$p, AOQ, type="l", 
ylab="AOQ", 
xlab="Fraction Nonconforming p")
title(paste("AOQL = ", formatC(max(AOQ))))
}
return(results)
}


print.SeqSampPlan=function(x,...){
print(data.frame(h1=x$h1, h2=x$h2, s=x$s))
}

plot.SeqSampPlan=function(x,y=NULL,...){
plot(x$k, x$accept, type="l", ylab=expression(d[k]), xlab="k", ylim=c(min(x$accept), max(x$reject)))
par(new=TRUE)
plot(x$k, x$reject, type="l", ylab="", xlab="", ylim=c(min(x$accept), max(x$reject)))
title("Sequential Acceptance Chart")
axis(1, tck = 1, col = "grey", lty = "dotted")
axis(2, tck = 1, col = "grey", lty = "dotted")
text(median(x$k), min(x$accept), "ACCEPT")
text(median(x$k), max(x$reject), "REJECT")
text(median(x$k), max(x$accept), "CONTINUE")
}
