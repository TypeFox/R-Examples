prediction_multiL_ff_norm <-
function(times,reg,ksb,ks1,ks21,ks22,kd1,kd2,delta1,delta2,noise)

{  

time2=c(-max(delta1,delta2),times)

h1=splinefun(time2,c(as.numeric(reg)[1],as.numeric(reg)))

p=c(ksb=ksb,ks1=ks1,ks21=ks21,ks22=ks22,kd1=kd1,kd2=kd2,delta1=delta1,delta2=delta2)



 sys<-function(t,z,p)

{ gp=p["ksb"]+p["ks1"]*h1(t-p["delta1"])-p["kd1"]*z

 list(gp,c())}



targ0=(p["ksb"]+p["ks1"]*h1(times[1]))/p["kd1"]

S1=deSolve::lsoda(targ0,times,sys,p,rtol=1e-4,atol=1e-6)

time2=c(-max(delta1,delta2),times)

reg2= S1[,2]

reg2= norm.data(reg2)+rnorm(n=length(times),mean=0,sd=abs(norm.data(reg2)*noise))



h2=splinefun(time2,c(as.numeric(reg2)[1],reg2))



sys<-function(t,z,p)

{ gp=p["ksb"]+p["ks21"]*h1(t-p["delta1"])+p["ks22"]*h2(t-p["delta2"])-p["kd2"]*z

list(gp,c())}



targ0=(p["ksb"]+p["ks21"]*h1(times[1])+p["ks22"]*h2(times[1]))/p["kd2"]

S2=deSolve::lsoda(targ0,times,sys,p,rtol=1e-4,atol=1e-6)



tar2 = S2[,2]

tar2 = norm.data(tar2)+rnorm(n=length(times),mean=0,sd=abs(norm.data(tar2)*noise))



output=list(tar1=reg2,tar2=tar2)

return(output)}
