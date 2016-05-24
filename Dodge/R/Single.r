
SSPlanBinomial=function(N,n,Ac, p=seq(0, 0.3, .001), Plots=TRUE)
{
OC = pbinom(Ac, n, p)
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

SSPlanHyper=function(N,n,Ac, p=seq(0, 0.3, .001), Plots=TRUE)
{
OC = phyper(Ac, N*p, N*(1-p), n)
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

SSPlanPoisson=function(N, n,Ac, p=seq(0, 0.3, .001), Plots=TRUE)
{
OC = ppois(Ac,n*p)
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

SSPDesignBinomial =function(AQL, alpha, LQL, beta){
nl=function(Ac, LQL, beta)
{
n=1
while(pbinom(Ac, n, LQL) >= beta){
n=n+1
} 
n
}
#nl(5, 0.04, 0.05)
#pbinom(5, 261, .04)

nu=function(Ac,AQL,alpha)
{
n=1
while(pbinom(Ac, n, AQL) >= 1-alpha){
n=n+1
} 
n
}
#nl(5, 0.01, 0.05)
#pbinom(5, 1049, .01)

Ac =0
while(nl(Ac, LQL, beta)>nu(Ac,AQL,alpha)){
Ac=Ac+1
}
n=nl(Ac, LQL, beta)
return(data.frame(n, Ac))
}

SSPDesignPoisson =function(AQL, alpha, LQL, beta)
{
nl=function(Ac, LQL, beta)
{
n=1
while(ppois(Ac, n* LQL) >= beta){
n=n+1
} 
n
}
#nl(5, 0.04, 0.05)
#ppois(5, 263* .04)

nu=function(Ac,AQL,alpha)
{
n=1
while(ppois(Ac, n* AQL) >= 1-alpha){
n=n+1
} 
n
}
#nl(5, 0.01, 0.05)
#ppois(5, 1052*.01)

Ac =0
while(nl(Ac, LQL, beta)>nu(Ac,AQL,alpha)){
Ac=Ac+1
}
n=nl(Ac, LQL, beta)
return(data.frame(n, Ac))
}

