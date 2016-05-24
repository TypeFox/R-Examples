prediction_doubleL_norm <-
function(times,reg1,reg2,ksb,ks1,ks2,kd,delta1,delta2,noise)

{   

   if (max(delta1,delta2)>0)

   {time2=c(-max(delta1,delta2),times)

   h1=splinefun(time2,c(as.numeric(reg1)[1],as.numeric(reg1)))

   h2=splinefun(time2,c(as.numeric(reg2)[1],as.numeric(reg2)))

   }else

   {h1=splinefun(times,as.numeric(reg1))

    h2=splinefun(times,as.numeric(reg2))}



  

         p=c(ksb=ksb,ks1=ks1,ks2=ks2,kd=kd,delta1=delta1,delta2=delta2)



    sys<-function(t,z,p)

   { gp=p["ksb"]+p["ks1"]*h1(t-p["delta1"])+p["ks2"]*h2(t-p["delta2"])-p["kd"]*z

     list(gp,c())

         }

   targ0=(p["ksb"]+p["ks1"]*h1(times[1])+p["ks2"]*h2(times[1]))/p["kd"]



   S=deSolve::lsoda(targ0,times,sys,p,rtol=1e-4,atol=1e-6)

   output=norm.data(S[,2])+rnorm(n=length(times),mean=0,sd=abs(norm.data(S[,2])*noise))



return(output)}
