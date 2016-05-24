powerfamily_twoparameter_balancevector=function(Deltas,requiredtypeIerror,requiredpower,delta0,delta1,sigma,J,weights)
{
     return(.C("powerfamily_twoparameter_nonintegern",as.double(Deltas[1]),as.double(Deltas[2]),as.double(requiredtypeIerror),as.double(requiredpower),as.double((delta1-delta0)/sigma),as.double(1),as.integer(J),finalparameters=double(3),ess=double(1),error=double(1),as.integer(4),as.double(weights))$ess)
}


extendedpowerfamily=function(Deltas,requiredtypeIerror,requiredpower,delta0,delta1,sigma,J)
{
design=.C("extendedpowerfamily",as.double(Deltas[1]),as.double(Deltas[2]),as.double(requiredtypeIerror),as.double(requiredpower),as.double((delta1-delta0)/sigma),as.double(1),as.integer(J),finalparameters=double(3))$finalparameters


temp=.C("operatingcharacteristics",as.double(Deltas[1]),as.double(Deltas[2]),as.double(design[2]),as.double(design[1]),logC2=double(1),as.double((delta1-delta0)/sigma),as.double(1),as.integer(J),typeIerror=double(1),power=double(1),finalparameters=double(2*J+1),ess_null=double(1),ess_crd=double(1),ess_dm=double(1))

return(list(groupsize=temp$finalparameters[1],futility=temp$finalparameters[seq(2,length(temp$finalparameters),by=2)],efficacy=temp$finalparameters[seq(3,length(temp$finalparameters),by=2)],ess=c(temp$ess_null,temp$ess_crd,temp$ess_dm),typeIerror=temp$typeIerror,power=temp$power,Deltas=Deltas,Ce=exp(design[2]),Cf=exp(temp$logC2)))



}


getparameters_twoparameter=function(Deltas,requiredtypeIerror,requiredpower,delta0,delta1,sigma,J,methodnumber,balanceparameter)
{
return(.C("powerfamily_twoparameter_nonintegern",as.double(Deltas[1]),as.double(Deltas[2]),as.double(requiredtypeIerror),as.double(requiredpower),as.double((delta1-delta0)/sigma),as.double(1),as.integer(J),finalparameters=double(3),ess=double(1),error=double(1),as.integer(methodnumber),as.double(balanceparameter))$finalparameters)
}

getparameters_twoparameter_balancevector=function(Deltas,requiredtypeIerror,requiredpower,delta0,delta1,sigma,J,weights)
{
return(.C("powerfamily_twoparameter_nonintegern",as.double(Deltas[1]),as.double(Deltas[2]),as.double(requiredtypeIerror),as.double(requiredpower),as.double((delta1-delta0)/sigma),as.double(1),as.integer(J),finalparameters=double(3),ess=double(1),error=double(1),as.integer(4),as.double(weights))$finalparameters)
}

powerfamily_fixedn_balancevector=function(par,requiredtypeIerror,requiredpower,requiredn,delta0,delta1,sigma,J,penaltyfactor,weights)
{
return(.C("powerfamily_fixedn",as.double(par[1]),as.double(par[2]),as.double(par[3]),as.double(requiredtypeIerror),as.double(requiredpower),as.double(requiredn),as.double((delta1-delta0)/sigma),as.double(1),as.integer(J),finalparameters=double(2*J+1),functionvalue=double(1),penaltyfactor=penaltyfactor,as.integer(4),as.double(weights))$functionvalue)
}

getoperatingcharacteristics=function(par,requiredn,delta0,delta1,sigma,J){
temp=.C("operatingcharacteristics",as.double(par[1]),as.double(par[2]),as.double(par[3]),as.double(requiredn),logC2=double(1),as.double((delta1-delta0)/sigma),as.double(1),as.integer(J),typeIerror=double(1),power=double(1),finalparameters=double(2*J+1),ess_null=double(1),ess_crd=double(1),ess_dm=double(1))

return(list(groupsize=temp$finalparameters[1],futility=temp$finalparameters[seq(2,length(temp$finalparameters),by=2)],efficacy=temp$finalparameters[seq(3,length(temp$finalparameters),by=2)],ess=c(temp$ess_null,temp$ess_crd,temp$ess_dm),typeIerror=temp$typeIerror,power=temp$power,Deltaf=par[1],Deltae=par[2],Cf=exp(par[3]),Ce=exp(temp$logC2)))

}

getexpectedsamplesizes=function(parameters,delta,sigma,J)
{
return(.C("getexpectedsamplesizes",as.double(parameters),as.double(delta),as.double(sigma),as.integer(J),expectedsamplesizes=double(50000))$expectedsamplesizes)
}

plotexpectedsamplesize=function(groupsize,futility,efficacy,delta,sigma,J,singlestagesamplesize,ylim,...)
{
   
parameters=rep(0,J*2+1)

parameters[1]=groupsize;
parameters[seq(2,length(parameters),by=2)]=futility
parameters[seq(3,length(parameters),by=2)]=efficacy

ess=getexpectedsamplesizes(parameters,delta,sigma,J)

#plot(seq(0,2*delta,length=50000),ess,type="l",xlab="Standardised treatment effect",ylab="Expected sample size per arm",lwd=2,col="blue",cex.axis=1.4,cex.lab=1.4,ylim=c(0,singlestagesamplesize))
#plot(x=seq(0,2*delta,length=50000),y=ess,type="l",xlab="Standardised treatment effect",ylab="Expected sample size per arm",lwd=2,col="blue",cex.axis=1.4,cex.lab=1.4,ylim=c(0,singlestagesamplesize))

 if(is.null(ylim))
     {ylim=c(0,1.1*singlestagesamplesize)}
plot(x=seq(0,2*delta,length=50000),y=ess,ylim=ylim,...)

abline(a=singlestagesamplesize,b=0,lty=2,lwd=2)

}


###############
###function to convert from known variance stopping boundaries to unknown 
###variance stopping boundaries
###############

converttounknownvarianceboundaries=function(design)
{

J=length(design$futility)
pvalues_futility=pnorm(design$futility)
pvalues_efficacy=pnorm(design$efficacy)
degreesoffreedom=(2*(1:J)*design$groupsize-2)

futility_unknownvariance=qt(pvalues_futility,degreesoffreedom)
efficacy_unknownvariance=qt(pvalues_efficacy,degreesoffreedom)

return(list(groupsize=design$groupsize,futility=futility_unknownvariance,efficacy=efficacy_unknownvariance,ess=design$ess,typeIerror=design$typeIerror,power=design$power,Deltaf=design$Deltaf,Deltae=design$Deltae,Cf=design$Cf,Ce=design$Ce))

}


###############
###Main function
###############




optgs=function(delta0=0,delta1=1/3,J=2,sigma=1,sd.known=TRUE,alpha=0.05,power=0.9,weights=c(0.95,0,0,0.05),initial=NULL)
{

#error checking:

#are delta0, delta1 and sigma real numbers with delta1>delta0 and sigma positive?
if(!is.double(delta0) || !is.double(delta1) || !is.double(sigma))
{
stop("delta0, delta1 and sigma must be real numbers")
}
if(delta0>=delta1)
{
stop("delta1 must be greater than delta0")
}

#is J an integer between 2 and 10?
if(ceiling(J)!=floor(J) || J<2 || J>10)
{
stop("Number of stages (J) must be an integer between 2 and 10")
}

#is weights of length 4 with all entries positive
if(length(weights)!=4 || !all(weights>=0) || sum(weights[1:3])==0)
{
stop("weights must be of length 4, all entries must be non-negative, and one of the first three entries must be strictly positive")
}
 
if(length(initial)==2)
{
if(initial[1]<(-0.5) | initial[1]>0.5 | initial[2]<(-0.5) | initial[2]>0.5)
{
stop("initial must be between -0.5 and 0.5")
}
}
if(length(initial)!=0 & length(initial)!=2)
{
stop("initial must be NULL or two-dimensional vector with entries between -0.5 and 0.5")
} 

if(!is.logical(sd.known))
{
stop("sd.known must be T or F")
}


weights=weights/sum(weights)

if(length(initial)==0)
{
initial=weights[1]*c(0.4,-0.4)+weights[2]*c(-0.2,0.4)+weights[3]*c(0.3,0.3)+weights[4]*c(-0.5,-0.5)
}

singlestagesamplesize=(2*sigma^2*(qnorm(1-alpha)+qnorm(power))^2)/(delta1-delta0)^2

optimaldeltas_nonintegern=optim(initial,powerfamily_twoparameter_balancevector,requiredtypeIerror=alpha,requiredpower=power,delta0=delta0,delta1=delta1,sigma=sigma,J=J,weights=weights,control=list(reltol=5e-5))


optimaldesign_nonintegern=getparameters_twoparameter_balancevector(optimaldeltas_nonintegern$par,requiredtypeIerror=alpha,requiredpower=power,delta0=delta0,delta1=delta1,sigma=sigma,J=J,weights=weights)


optimaln_floor=floor(optimaldesign_nonintegern[1])

optimal_floor=optim(c(optimaldeltas_nonintegern$par,optimaldesign_nonintegern[2]),powerfamily_fixedn_balancevector,requiredtypeIerror=alpha,requiredpower=power,requiredn=optimaln_floor,delta0=delta0,delta1=delta1,sigma=sigma,J=J,penaltyfactor=100*singlestagesamplesize,weights=weights,control=list(reltol=1e-8))

optimaln_ceil=ceiling(optimaldesign_nonintegern[1])

optimal_ceil=optim(c(optimaldeltas_nonintegern$par,optimaldesign_nonintegern[2]),powerfamily_fixedn_balancevector,requiredtypeIerror=alpha,requiredpower=power,requiredn=optimaln_ceil,delta0=delta0,delta1=delta1,sigma=sigma,J=J,penaltyfactor=100*singlestagesamplesize,weights=weights,control=list(reltol=1e-8))

if(optimal_ceil$value>optimal_floor$value)
{
optimaldesign=optimal_floor;
optimaln=optimaln_floor
}

if(optimal_ceil$value<=optimal_floor$value)
{
optimaldesign=optimal_ceil;
optimaln=optimaln_ceil
}


trialproperties=getoperatingcharacteristics(optimaldesign$par,optimaln,delta0,delta1,sigma,J)

if(!sd.known)
{
trialproperties=converttounknownvarianceboundaries(trialproperties)
}

trialproperties$J=J
trialproperties$delta0=delta0
trialproperties$delta1=delta1
trialproperties$sigma=sigma
trialproperties$singlestageSS=singlestagesamplesize

class(trialproperties)="OptGS"

return(trialproperties)


}


###############
###Generic print function
###############


print.OptGS=function (x,...)
    {

cat("\nGroupsize: ", paste(deparse(round(x$groupsize,digits=1)), sep = "\n", collapse = "\n"),"\nFutility boundaries ",paste(round(x$futility,digits=2),sep=" "),"\nEfficacy boundaries ",paste(round(x$efficacy,digits=2),sep=" "),"\nESS at null:    ",paste(round(x$ess[1],digits=1)),sep=" ","\nESS at CRD:     ",paste(round(x$ess[2],digits=1),sep=" "),"\nMaximum ESS:    ",paste(round(x$ess[3],digits=1),sep=" "),"\nMax sample-size:",paste(round(x$J*x$groupsize,digits=1),sep=" "),"\n")
}


###############
###Prints expected sample size of the required design
##############

plot.OptGS=function(x,ylim=NULL,...)
{
plotexpectedsamplesize(x$groupsize,x$futility,x$efficacy,x$delta1-x$delta0,x$sigma,x$J,x$singlestageSS,ylim,...)
}



##############
###Returns extended power-family design for user-specified (Deltaf,Deltae)
##############

powerfamily=function(futility=0,efficacy=0,delta0=0,delta1=1/3,J=2,sigma=1,sd.known=TRUE,alpha=0.05,power=0.9)
{
if(!is.double(delta0) || !is.double(delta1) || !is.double(sigma))
{
stop("delta0, delta1 and sigma must be real numbers")
}
if(delta0>=delta1)
{
stop("delta1 must be greater than delta0")
}

#is J an integer between 2 and 10?
if(ceiling(J)!=floor(J) || J<2 || J>10)
{
stop("Number of stages (J) must be an integer between 2 and 10")
}


if(!is.logical(sd.known))
{
stop("sd.known must be T or F")
}

if(!is.double(efficacy) || efficacy<= -1 || efficacy>1)
{
stop("Efficacy shape parameter must be between -1 and 1")
}

if(!is.double(futility) || futility<= -1 || futility>1)
{
stop("Futility shape parameter must be between -1 and 1")
}

singlestagesamplesize=(2*sigma^2*(qnorm(1-alpha)+qnorm(power))^2)/(delta1-delta0)^2



trialproperties=extendedpowerfamily(Deltas=c(futility,efficacy),alpha,power,delta0,delta1,sigma,J)


if(!sd.known)
{
trialproperties=converttounknownvarianceboundaries(trialproperties)
}

trialproperties$J=J
trialproperties$delta0=delta0
trialproperties$delta1=delta1
trialproperties$sigma=sigma
trialproperties$singlestageSS=singlestagesamplesize
trialproperties$Deltaf=futility
trialproperties$Deltae=efficacy
trialproperties$Cf=trialproperties$Cf
trialproperties$Ce=trialproperties$Ce



class(trialproperties)="OptGS"



return(trialproperties)

}
