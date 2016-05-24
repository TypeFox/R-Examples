#
# vim:set ff=unix expandtab ts=2 sw=2:
test.MeanAge3=function(){
   pf=function(str){print(paste(str,"=",eval(parse(text=str))))}
   # we consider a model with input and outputrates as functions of time
   # described by a possibly nonlinear ode.
   # First we compute the solution of the system.

   C0=c(5,150)
   k1=1/100
   k2=1/300
   GE=read.csv("GlobalEmissions.csv")[c("Year","FossilFuel")]
   YearsMeasured=GE[,1]
   firstYear=min(YearsMeasured)
   lastYear=max(YearsMeasured)
   FF=GE[,2]
   firstFF=FF[1]
   lastFF=FF[length(FF)]
   minFF=min(FF)
   spFF=splinefun(YearsMeasured,FF)
   firstApprox=function(t){
      perturbation=minFF*(1+0.5*sin(2*pi/4*t))
      if (t<firstYear){CInTeraGramm=firstFF+perturbation}
      #if (t>lastYear){CInTeraGramm=firstFF}
      if (t>lastYear){CInTeraGramm=lastFF}
      if (t>=firstYear & t<=lastYear){CInTeraGramm=spFF(t)}
      # the dataset contains the mass of C (not CO2) in the atmosphere in units of Teragramms
      # to include it in the model we have to translate it to molar numbers
      # which we do by means of the molar weight of Carbon =12.011g/mol

      CInmols=CInTeraGramm*1e1/12.011
      # to avoid numerical instability triggered by the loss of significance in the computaion of the linear operator
      # for steady state solutions (where the derivative become 0) we perturb the system a litle bit in the initial phase

      return(CInmols)
   }
   times=seq(1,2300)
   unsmoothed=mapply(firstApprox,times)
   model=smooth.spline(times,unsmoothed,spar=0.0)
   CO2inputrate=function(t){predict(model,t)$y+200}
   pdf(file="runit.MeanAge3.pdf",paper="a4r")
   c=c("black","red","green","blue")
   lts=c(1,2)
   lws=c(8,4)
   plot(times,CO2inputrate(times),type="l",lwd=lws[1],col=c[1],lty=lts[1])
   I0=k1*C0[1]
   Idot1=function(Y,t){CO2inputrate(t)}
   #Idot1=function(Y,t){I0*(1+0.9*sin(t/50))}
   Idot2=function(Y,t){-Odot1(Y,t)/2}
   Idot=function(Y,t){matrix(byrow=TRUE,c(Idot1(Y,t),Idot2(Y,t)))}
   
   Odot1=function(Y,t){-k1*Y[1]^2}
   Odot2=function(Y,t){-k2*Y[1]*Y[2]}
   Odot=function(Y,t){matrix(byrow=TRUE,c(Odot1(Y,t),Odot2(Y,t)))}
   
   Ydot=function(Y,t){Idot(Y,t)+Odot(Y,t)}
   tstart=0
   tend=210
   tn=1000
   tol=.02/tn
   maxage=tend-tstart
   times=seq(tstart,tend,maxage/tn)
   sol=solver(times,Ydot,C0)
   fs1=splinefun(times,sol[,1])
   fs2=splinefun(times,sol[,2])
   fs=function(t){matrix(nrow=2,
	mapply(function(fun){fun(t)},list(fs1,fs2))
   )}
   ######################################################################
   # With the help of the solution we can re express 
   # the system as a linear problem with the same solution 
   # which enables us to track normalized amounts
   # of matter trough the system.
   OdotLin1=linMaker(Odot1,fs,fs1)
   OdotLin2=linMaker(Odot2,fs,fs2)
   
   IdotLin1=linMaker(Idot1,fs,fs1)
   IdotLin2=linMaker(Idot2,fs,fs2)
   
   IdotT1=function(t){Idot1(fs(t),t)}
   IdotT2=function(t){Idot2(fs(t),t)}
   
   CdotLin1=function(Y,t){IdotLin1(Y,t)+OdotLin1(Y,t)}
   CdotLin2=function(Y,t){IdotLin2(Y,t)+OdotLin2(Y,t)}
   fsLin1=splinefun(
      times,
      solver( times, CdotLin1, C0[1])
   )
   fsLin2=splinefun(
      times,
      solver( times, CdotLin2, C0[2])
      # solver( times, IdotLin2,0)
      #+solver( times, OdotLin2,C0[2])
      )
   plot(times,fs1(times),type="l",lwd=lws[1],col=c[1],lty=lts[1])
   lines(times,fsLin1(times),col=c[2],lty=lts[2],lwd=lws[2])
   legend(
   "bottomleft",
     c( "solution ","solution of the equivalent linear system"),
     lty=lts,
     col=c(c[1],c[2])
   )
   plot(times,fs2(times),type="l",lwd=lws[1],col=c[1],lty=lts[1])
   lines(times,fsLin2(times),lwd=lws[2],col=c[2],lty=lts[2])
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
   #ma3=splinefun(times,MeanAge3(IdotT,OdotLin,fs,times))
   #plot(times,ma3(times))
   dev.off()
}
