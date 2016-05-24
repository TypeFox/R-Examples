#####################################################################################
####################################################################################
#
#	profile likelihood for etasclass objects
#
#
#####################################################################################
#####################################################################################
profile.etasclass<-function(fitted,iprofile		=4,
				 nprofile	=7,
				 kprofile	=3,
				 profile.approx	=FALSE,...){
if(class(fitted)!="etasclass") stop("object is not of the required class etasclass")
iprofile=trunc(iprofile)
if(iprofile>8 | iprofile <1 ) stop("'iprofile' must be an integer between 1 and 8")

nprofile=trunc(nprofile)
if(nprofile<1) stop("'nprofile' must be at least 1")
fitted$params.ind=as.numeric(fitted$params.ind)
if(fitted$params.ind[iprofile]==0) stop("cannot compute profile for a parameter not estimated")


		logl.vec	=0
		param.vec	=0
		
		params.indprof	=fitted$params.ind
		params.fix	=fitted$params.fix
		n.params	=sum(params.indprof)
		iterlim		=100
		n		=length(fitted$cat$time)
#####################################################################################
if(profile.approx){
#	approximation of second order for initial values for non-profile estimators			

		ind.prof	=sum(params.indprof[1:iprofile])
		ind.noprof	=setdiff(1:n.params,ind.prof)
		delta.psi	=solve(fitted$risult$hessian[ind.noprof,ind.noprof])%*%fitted$risult$hessian[ind.noprof,ind.prof]
}
#####################################################################################

		params.indprof[iprofile]	=0
		sq			=fitted$sqm[iprofile]*kprofile
		
		param.vec		=seq(fitted$params[iprofile]-sq,fitted$params[iprofile]+sq,length.out=nprofile)
		logl.vec		=array(0,nprofile)
		params.MLtot		=fitted$params
		
		for (j in 1:nprofile)	{

	cat("start profile-L computation. j=  ")
  cat(j,"\n")

		params.fix[iprofile]	=param.vec[j]
		params		=log(params.MLtot[params.indprof==1])
	
		if(profile.approx) params=params-delta.psi*(params.fix[iprofile]-params.MLtot[iprofile])
cat(params,"\n")
if (fitted$usenlm){
	risult.profile =nlm(etas.mod2,params,
#		hessian	=TRUE,
#		stepmax=200,
#		steptol	=1e-6,
		typsize=abs(params),
		#fscale=l.optim,
		ndigit=5,print.level=0,
		iterlim		=100,
		params.ind=params.indprof,
 		params.fix=params.fix,
		cat=fitted$cat,
		magn.threshold=fitted$magn.threshold,
		back.dens=fitted$back.dens,
		back.integral=fitted$back.integral,
		onlytime=fitted$onlytime,
		rho.s2=fitted$rho.s2,
		ntheta=fitted$ntheta)

	}
else {
	risult.profile =optim(params,etas.mod2,
#		stepmax=200,
#		steptol	=1e-6,
		method		="BFGS",
		hessian		=TRUE,
		control		=list(maxit=iterlim,fnscale=n/diff(range(cat$time)),parscale=sqrt(exp(params))),
		params.ind=params.indprof,
 		params.fix=params.fix,
		cat=fitted$cat,
		magn.threshold=fitted$magn.threshold,
		back.dens=fitted$back.dens,back.integral=fitted$back.integral,
		onlytime=fitted$onlytime,
		rho.s2=fitted$rho.s2,
		ntheta=fitted$ntheta)
	}

logl.vec[j]=ifelse(fitted$usenlm,risult.profile$minimum,risult.profile$value)
cat(" profile likelihood found","\n")
cat(c(j,param.vec[j],logl.vec[j]),"\n")
	}
			
			ris=list(
			iprofile	=iprofile,
			logl.vec	=logl.vec,
			param.vec	=param.vec,
			general.optimum =c(fitted$params[iprofile],fitted$logl)
			)

			
class(ris)=c("profile.etasclass",class(ris))

return(ris)


#############################################################################################
# end of profile computation
#############################################################################################

#############################################################################################
}
