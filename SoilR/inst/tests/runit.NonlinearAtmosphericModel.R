#
# vim:set ff=unix expandtab ts=2 sw=2:
testSolution=function(){
  require(parallel)
   k12=0.069565
   k21=0.0136298e-19
   k23=0.055555
   k31=0.0095236e-19

   #F1=3.75575e19 #total O2=CO2+O2
   #F2=153.333e15 #total C=C2+C3+CO2
   
   # To create the inputrate function we use real fossil fuel emission data for 
   # the time where they are available and extrapolate them as >> constant << before
   # and after the times measured
   GE=read.csv("GlobalEmissions.csv")[c("Year","FossilFuel")]
   YearsMeasured=GE[,1]
   firstYear=min(YearsMeasured)
   lastYear=max(YearsMeasured)
   FF=GE[,2]
   firstFF=FF[1]
   lastFF=FF[length(FF)]
   minFF=min(FF)
   pp("minFF",environment())
   spFF=splinefun(YearsMeasured,FF)
   convfact= 1e12/12.011
   firstApprox=function(times){
      h=function(t){
	 perturbation=minFF*(1+0.9*cos(2*pi/200*t))
      	 if (t<firstYear){CInTeraGramm=firstFF+perturbation}
      	 if (t>lastYear){CInTeraGramm=firstFF+perturbation}
      	 #if (t>lastYear){CInTeraGramm=lastFF}
      	 if (t>=firstYear & t<=lastYear){CInTeraGramm=spFF(t)+perturbation}
      	 # the dataset contains the mass of C (not CO2) in the atmosphere in units of Teragramms
      	 # to include it in the model we have to translate it to molar numbers
      	 # which we do by means of the molar weight of Carbon =12.011g/mol

	 CInmols=CInTeraGramm*convfact
	 return(CInmols)
      }
      return(mapply(h,times))
   }
   times=seq(0,2300,00.1)
   unsmoothed=firstApprox(times)
   model=smooth.spline(times,unsmoothed,spar=0.3)
   smoothed=predict(model,times)$y
   firstApproxSmoothed=function(t){predict(model,t)$y}
   CO2inputrate=firstApprox
   #CO2inputrate=firstApproxSmoothed

   Idot4=function(Y,t){
      C2 =Y[3]
      return(k23*C2)
   }
   # we define the outputs of the 4 reservoirs as functions of 
   Odots=c(
	function(Y,t){-k12*Y[1]},
	function(Y,t){-k31*Y[4]*Y[2]-k21*Y[3]*Y[2]},
	function(Y,t){-k21*Y[3]*Y[2]-k23*Y[3]},
	function(Y,t){-k31*Y[4]*Y[2]}
	)
   Odot=vecFuncMaker(Odots,Y,t)
   Idots=c(
 	function(Y,t){k31*Y[4]*Y[2]+k21*Y[3]*Y[2]+CO2inputrate(t)},
 	function(Y,t){k12*Y[1]+startValues[2]*0*(cos(2*pi/200*t+pi/4))},
 	function(Y,t){k12*Y[1]+startValues[3]*0*(cos(2*pi/300*t-pi/4))},
 	function(Y,t){k23*Y[3]+startValues[4]*0*(cos(2*pi/100*t+pi/3))}
      )
   Idot=vecFuncMaker(Idots,Y,t)
   ydot=function(Y,t){
      Idot(Y,t)+Odot(Y,t)
   }
   tstart=0
   tend=2100
   maxage=tend-tstart
   tn=2100
   increment=maxage/tn
   #times=seq(tstart,tend)
   times=seq(tstart,tend,increment)
   startValues=c(
      57.5e15,
      3.75e19,
      37.5e15,
      58.333e15)
   sol=solver(times,ydot,startValues)
   ss=mclapply(
	#create a list of lists where every sublist contains a column of the matrix
	mcmapply(function(i){list(sol[,i])},1:ncol(sol)), 
	#create a spline aproximation for every column
	function(vec){splinefun(times,vec)}
   )	
   sY=function(t){matrix(nrow=length(ss),
	mapply(function(fun){fun(t)},ss)
   )}
   # get new startvalues somewhere away from the steady state
   tm=tstart
   #tm=1700
   newStartValues=sY(tm)
   lintimes=seq(tm,tend,increment)

   # create the linearized versions of the operator parts
   # outputs
   OdotLins=mapply(function(i){linMaker(Odots[[i]],sY,ss[[i]])},1:length(Odots))
   # inputs
   IdotLins=mapply(function(i){linMaker(Idots[[i]],sY,ss[[i]])},1:length(Idots))
   # complete operator
   CdotLins=mapply(
        function(i){
        	function(Y,t){IdotLins[[i]](Y,t)+OdotLins[[i]](Y,t)}
        },	
        1:length(OdotLins)
   )
   checksols=mapply(
	function(i){
	   splinefun(
	       lintimes,
	       solver(lintimes,CdotLins[[i]],newStartValues[i])
	   )
	},
	1:length(CdotLins)
   )
   Idot4T=splinefun(times,sapply(times,function(t){Idot4(sY(t),t)}))
   
   
   #create plots
   pdf(file="runit.NonlinearAtmosphericModelLong.pdf",paper="a4r")
   ptimes=seq(1600,1800) 
   #plot(ptimes,CO2inputrate(ptimes),type="l",lty=1,lwd=4,col="red")
   #lines(ptimes,firstApproxSmoothed(ptimes),type="l",lty=2,lwd=4,col="black")
  #mapply(
  #       function(i){ 
  #          plot(lintimes,ss[[i]](lintimes),type="l",col="black",lty=1)
  #          lines(lintimes,checksols[[i]](lintimes),col="red",lty=2)
  #       },
  #       seq(1,4)
  # )

   meanAge=splinefun(times,MeanAge2(Idot4T,OdotLins[[4]],checksols[[4]],times))
  # meanTT=splinefun(MeanTT(OdotLins[[4]],times))
   lt1=1
   c=seq(1,5)
   par(mar=c(5,4,4,5)+.1)
   ptimes=seq(1800,2100)
   CO2=sol[,1]
   O2=sol[,2]
   C2=sol[,3]
   C3=sol[,4]
   plot (ptimes,ss[[1]](ptimes),type="l",lwd=4,lty=lt1,col=c[1],
        ylab="amount of Carbon in mol",
        xlab="Time",
        ylim=c(
        	#min(CO2,C2,C3),
        	0,
        	max(CO2,C2,C3)
        ))
   #lines(times,ss[[2]](ptimes),type="l",lty=lt1,col=c[2])
   lines(ptimes,ss[[3]](ptimes),type="l",lwd=4,lty=lt1,col=c[3])
   lines(ptimes,ss[[4]](ptimes),type="l",lwd=4,lty=lt1,col=c[4])
   #lines(ptimes,checksols[[1]](ptimes),type="l",lty=2,col="red",lwd=4)
   #lines(ptimes,checksols[[3]](ptimes),type="l",lty=2,col="green",lwd=4)
   lines(ptimes,checksols[[4]](ptimes),type="l",lty=2,col="yellow",lwd=4)
   par(new=TRUE)
   plot (ptimes,sapply(ptimes,CO2inputrate),type="l",lty=lt1,lwd=4,col=c[2],xaxt="n",yaxt="n",xlab="",ylab="")
   axis(4)
   mtext("Carbon input by fossil fuel combustion in mol/year",side=4,line=3)
   legend(
   "topleft",
     c( "CO2","C2","C3","FF"),
     lty=rep(lt1,4),
     col=c(c[1],c[3],c[4],c[2])
   )
  # plot (ptimes,meanTT(ptimes),type="l",lty=1,col=c[1],ylab=expression("mean transit times"),xlab="Time")
  # lines(ptimes,meanAge(ptimes),type="l",lty=2,col=c[2],ylab=expression("mean age"),xlab="Time")
  #plot (times,CO2,type="l",lty=lt1,col=c[1],ylab=expression("CO[2] in 1e15 mol"),xlab="Time",ylim=c(0,max(CO2)))
  #plot (times,O2,type="l",lty=lt1,col=c[2],ylab="O2 in 10^15 mol",xlab="Time",ylim=c(min(O2),max(O2)))
  #plot (times,C2,type="l",lty=lt1,col=c[3],ylab="C2 in 10^15 mol",xlab="Time",ylim=c(0,max(C2)))
  #plot (times,C3,type="l",lty=lt1,col=c[4],ylab="C3 in 10^15 mol",xlab="Time",ylim=c(0,max(C3)))
  #plot (times,sapply(times,CO2inputrate),type="l",lty=lt1,col=c[3],ylab="ExternalIdot",xlab="Time")
   dev.off()
   pdf(file="runit.NonlinearAtmosphericModelMeanAge.pdf",paper="a4r")
   plot(ptimes,meanAge(ptimes),type="l",lty=2,col=c[2],ylab=expression("mean age"),xlab="Time")
   dev.off()
}

