prediction_L_norm <-
function(times,reg,ksb,ks,kd,delta,noise)

{  

if (delta>0)

{time2=c(-delta,times)

h=splinefun(time2,c(as.numeric(reg)[1],as.numeric(reg)))

}else

{h=splinefun(times,reg)}



  

p=c(ksb=ksb,ks=ks,kd=kd,delta=delta)



sys<-function(t,z,p)

{ gp=p["ksb"]+p["ks"]*h(t-p["delta"])-p["kd"]*z

list(gp,c())}



targ0=(p["ksb"]+p["ks"]*h(times[1]))/p["kd"]



S=deSolve::lsoda(targ0,times,sys,p,rtol=1e-4,atol=1e-6)

tar = S[,2]

tar = norm.data(tar) +rnorm(n=length(times),mean=0,sd=abs(norm.data(tar)*noise))



return(tar)}
