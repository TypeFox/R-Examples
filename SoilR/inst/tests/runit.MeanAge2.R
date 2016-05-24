#
# vim:set ff=unix expandtab ts=2 sw=2:
test.MeanAge2=function(){
   pf=function(str){print(paste(str,"=",eval(parse(text=str))))}
   # we consider a model with input and outputrates as functions of time
   # described by a possibly nonlinear ode.
   # First we compute the solution of the system.

   C0=c(1/2)
   k=1/10
   I0=k*C0
   Odot=function(Y,t){-k*Y^2}
   #Odot=function(Y,t){-k*Y}
   Idot=function(Y,t){I0*(1+0.9*sin(t/50))}
   #Idot=function(Y,t){I0}
   Ydot=function(Y,t){Idot(Y,t)+Odot(Y,t)}
   tstart=0
   tend=2000
   tn=2000
   tol=.02/tn
   maxage=tend-tstart
   times=seq(tstart,tend,maxage/tn)
   fs=splinefun(times,solver(times,Ydot,C0))
   pdf(file="runit.MeanAge2.pdf",paper="a4r")
   c=c("black","red","green","blue")
   lts=c(1,2)
   lws=c(8,4)
   plot(times,fs(times),type="l",lwd=lws[1],col=c[1],lty=lts[1],ylim=c(0,max(fs(times))))
   ######################################################################
   # With the help of the solution we can re express 
   # the system as a linear problem with the same solution 
   # which enables us to track normalized amounts
   # of matter trough the system.
   OdotLin=function(Y,t){Y*Odot(fs(t),t)/fs(t)}
   IdotT=function(t){Idot(fs(t),t)}
   fsLin=splinefun(
      times,
      solver(
             times,
             function(Y,t){Idot(Y,t)+OdotLin(Y,t)},
             C0
      )
   )
   lines(times,fsLin(times),col=c[2],lty=lts[2],lwd=lws[2])
   legend(
   "bottomleft",
     c( "solution ","solution of the equivalent linear system"),
     lty=lts,
     col=c(c[1],c[2])
   )

   #ma=splinefun(times,MeanAge(IdotT,OdotLin,fs,times))
   #plot(times,ma(times))
   #ma2=splinefun(times,MeanAge2(IdotT,OdotLin,fs,times))
   #plot(times,ma2(times))
   ma3=splinefun(times,MeanAge3(IdotT,OdotLin,fs,times))
   plot(times,ma3(times))
   dev.off()
}
