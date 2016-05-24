#
# vim:set ff=unix expandtab ts=2 sw=2:
test.meanAge=function(){
   pf=function(str){print(paste(str,"=",eval(parse(text=str))))}
   # we consider a model with input and outputrates as functions of time
   # described by a possibly nonlinear ode.
   # The solution of the system has to be computed first and the input and output terms .

   C0=c(10)
   #model the startvalues by means of an inpulsive Input
   impulse=function(t){2*C0*dnorm(t,mean=0,sd=1e-8)}
   k=1/10
   I0=k*C0
   #Odot=function(Y,t){-k*Y^8}
   Odot=function(Y,t){-k*Y}
   #Idot=function(Y,t){I0*(1+0.9*sin(t/50))}
   Idot=function(Y,t){I0}#+impulse(t)}
   Ydot=function(Y,t){Idot(Y,t)+Odot(Y,t)}
   t0=0
   tend=100
   tn=10
   tol=.02/tn
   times=seq(t0,tend)
   fs=splinefun(times,solver(times,Ydot,C0))
   pdf(file="runit.MeanAge.pdf",paper="a4r")
   ######################################################################
   # With the help of the solution we can re express 
   # the system as a linear problem with the same solution 
   # which enables us to track normalized amounts
   # of matter trough the system.
   koeff=function(t){Odot(fs(t),t)/fs(t)}
   #plot(times,mapply(koeff,times),type="l",lwd=4,col="black",lty=1)
   OdotLin=function(Y,t){Y*Odot(fs(t),t)/fs(t)}
   IdotLin=function(Y,t){Y*Idot(fs(t),t)/fs(t)}
   YdotLin=function(Y,t){IdotLin(Y,t)+OdotLin(Y,t)}
   IdotT=function(t){Idot(fs(t),t)}
   fsLin=splinefun(
      times,
      solver(
             times,
             function(Y,t){Idot(Y,t)+OdotLin(Y,t)},
             C0
      )
   )
   ######################################################################
   anasol=C0 #steady state solution
   #plot(times,fs(times),type="l",lwd=4,col="black",lty=1)
   #lines(times,fsLin(times),type="l",lwd=4,col="red",lty=2)
   

   ##rho=function(a,t){
   ##t_i=t-a   
   ##sa=solver(t_i,t,YdotLin,IdotLin(t-a)
   ##   val=sa(t)/fs(t)
   ##   return(val)
   ##}
   ##ma2=splinefun(times,MeanAge2(IdotT,OdotLin,fs,times))

   ma3=splinefun(times,MeanAge3(IdotT,OdotLin,fs,times))
   ma4=splinefun(times,MeanAge3(function(t){0},OdotLin,fs,times))
   anaMean=function(t){
      return(t0*exp(-k*t0) - t*exp(-k*t) + (exp(-k*t0)-exp(-k*t))/k+exp(-t*k)*t)
   }
   anaMean2=function(t){
      return(t)
   }
   plot(times,ma3(times),type="l",lwd=4,col="black",lty=1,ylab="Mean Age",xlab="Time")
   lines(times,mapply(anaMean,times),type="l",lwd=4,col="red",lty=2)#,ylim=c(0,max(fs(times))))
   legend(
   "bottomright",
     c( "numerical solution ","anlytical solution"),
     lty=c(1,2),
     col=c("black","red")
   )
   dev.off()
   pdf(file="runit.MeanAgeImpuls.pdf",paper="a4r")
   plot(times,ma4(times),type="l",lwd=4,col="black",lty=1,ylab="Mean Age",xlab="Time")
   lines(times,mapply(anaMean2,times),type="l",lwd=4,col="red",lty=2)#,ylim=c(0,max(fs(times))))
   legend(
   "bottomright",
     c( "numerical solution ","anlytical solution"),
     lty=c(1,2),
     col=c("black","red")
   )

   dev.off()
}
