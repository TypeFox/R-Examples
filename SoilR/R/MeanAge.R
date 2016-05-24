#
# vim:set ff=unix expandtab ts=2 sw=2:
MeanAge=function# mean age for a general one pool model
### The function computes the mean age for one pool of a possibly nonlinear model
(IdotT,	##<< The inputrate as a function of time 
OdotLin,##<< The outputrate of this pool as a linear operator (a function of y and t)
sol,	##<< The solution of the nonlinear equation for this pool as a function of time
times 	##<< A vector containing the points in time where the solution is sought.
){
rho=function(a,tk){
	startTime=tk-a
	#startVal=IdotT(startTime)
	startVal=sol(startTime)
	if (startTime<0){
	        print("if")
		print(paste("a=",a,"tk=",tk))
		return(startVal)
	}
	times=c(startTime,tk,2*tk)
	s=solver(times,OdotLin,startVal)
	res=s[2]/sol(tk)
	return(res)
}
E=function(tk){
	E_integrand=function(Y,a){a*rho(a,tk)}
	res=solver(c(0,tk*0.99),E_integrand,0)[2]
return(res)
}

res=sapply(times,E)
return(res)
   ### A vector containing the mean age for the specified times
}
########################################################
MeanAge2=function# mean age for a general one pool model
### The function computes the mean age for one pool of a possibly nonlinear model
(IdotT,	##<< The inputrate as a function of time 
OdotLin,##<< The outputrate of this pool as a linear operator (a function of y and t)
sol,	##<< The solution of the nonlinear equation for this pool as a function of time
times 	##<< A vector containing the points in time where the solution is sought.
){
   maxage=max(times)-min(times)
fineTimes=(seq(sqrt(min(times)),sqrt(max(times)),sqrt(maxage)/10000))^2   
#fineTimes=seq(min(times),max(times),maxage/100000)   
startval=1
OofT=splinefun(fineTimes,solver(fineTimes,OdotLin,startval))
rho=function(a,tk){
	s=OofT(tk)*IdotT(tk-a)/OofT(tk-a)
	res=s/sol(tk)
	return(res)
}
E=function(tk){
	rho_tk=function(a){rho(a,tk)}
	E_integrand=function(Y,a){a*rho(a,tk)}
	res=solver(c(tk/1000,tk),E_integrand,0)[2]
return(res)
}

res=mcmapply(E,times,mc.cores=16)
return(res)
   ### A vector containing the mean age for the specified times
}
########################################################
MeanAge3=function# mean age for a general one pool model
### The function computes the mean age for one pool of a possibly nonlinear model
(IdotT,	##<< The inputrate as a function of time 
OdotLin,##<< The outputrate of this pool as a linear operator (a function of y and t)
sol,	##<< The solution of the nonlinear equation for this pool as a function of time
times 	##<< A vector containing the points in time where the solution is sought.
){
so=solver(times,OdotLin,1)
#we assume that the startvalue is the result of an input in the first timestep
# first we compute the value of the solution at the first time
# this is also very useful if the initial value is zero
# which would make our Linear Operator useless since it has a fixpoint at 0 (being linear)

deltaT=times[[2]]-times[[1]]
tstart=times[[2]] 
startvalue=sol(tstart)
if(startvalue==0){stop("cannot compute the solution since we are in a fixed point")}
# now we ad an impulsive inputrate necessarry to achieve the startvalue 
# nearly immidiately but possibly smooth for the integration

Id=function(t){IdotT(t)}
l=100

pdf(file="runit.MeanAge.Densities.pdf",paper="a4r")
#plot(times,so)
OofT=splinefun(times,so)
#plot(times,OofT(times))
s=function(a,tk){
   OofT(tk)*Id(tk-a)/OofT(tk-a)
}
rho=function(a,tk){
	res=s(a,tk)/sol(tk)
	return(res)
}
E=function(tk){
        l=100
	#ages=log(seq(1,exp(tk*(1-1/l)),length=l))
        ages=seq(0,tk-times[[2]],length=l)
	s_tk=function(a){s(a,tk)}
	ss_tk=splinefun(ages,mapply(s_tk,ages))
	# we compute the probability for a partikle to 
	# have an age in the range of [0,tk-times[[2]]]
	one=integrate(ss_tk,lower=0,upper=tk-times[[2]])$value/sol(tk)
	# so error=1-one must be the probability to have an age in the range
        #[tk-times[[2]],tk] so we estimate rho for this inteval as
	rho_0=(1-one)/(deltaT)
        Integrand=function(a){a*rho(a,tk)} 
        spIntegrand=splinefun(ages,sapply(ages,Integrand))
	res=integrate(spIntegrand,lower=0,upper=max(ages))$value+rho_0*(tk-(deltaT)/2)
	return(res)
}
times[[1]]<-times[[2]]
#times[[2]]<-times[[3]]
#print(times)
dev.off()
#lapply(times,E)
res=mclapply(times,E,mc.cores=16)
return(res)
   ### A vector containing the mean age for the specified times
}
