## package for fitting dependent mixture models

# .onLoad <- function(lib, pkg) {
# 	library.dynam("depmix", pkg, lib)
# }
# 
# .onUnLoad <- function(libpath) {
# 	library.dynam.unload("depmix",libpath)
# }

###################################
#                                 #
#   FIT DEPENDENT MIXTURE MODEL   #
#                                 #
###################################

# if(getRversion() >= "2.15.1")  utils::globalVariables(c("donlp2", "donlp2Control"))

fitdmm <- function(dat,dmm,printlevel=1,poster=TRUE,tdcov=0,ses=TRUE,method="optim",vfactor=15,der=1,iterlim=100,kmst=!dmm$st,kmrep=5,postst=FALSE) { 
	mod=dmm #keep copy of the original model
## create good starting values if not provdided
	if(kmst) {
		if(class(dmm)[1]!="dmm") {
			warning("No automatic starting values for other than single component 1-group models.")
			kmst=FALSE
		}
		if(nitems(dat)==1 && dmm$itemtypes[1]>1) {
			warning("No automatic starting values for univariate multinomial time series.")
			kmst=FALSE
		}
	}
	
	if(kmst) {
		if(dmm$nstates>1) {
			logl=matrix(0,kmrep,2)
			stvalskm=matrix(0,kmrep,dmm$npars)
			stvalspost=matrix(0,kmrep,dmm$npars)
			for(i in 1:kmrep) {
				stvalskm[i,]=kmstart(dat,dmm)
				dmm$pars=stvalskm[i,]
				logl[i,1]=loglike(dat,dmm)$logl
				if(postst) {
					stvalspost[i,]=poststart(dat,dmm)
					dmm$pars=stvalspost[i,]
					logl[i,2]=loglike(dat,dmm)$logl
				} else {
					stvalspost[i,]=stvalskm[i,]
					logl[i,2]=logl[i,1]
				}
				if(printlevel>9) cat("logl st model ", i, ": ", logl[i,1], "\n")
			}
			dmm$pars=stvalspost[which.max(logl[,2]),]
			dmm$st=TRUE
		}
	} # end else
	
	if(postst) {
		if(class(dmm)[1]!="dmm") {
			warning("No automatic boosting of starting values for other than single component 1-group models.")
			postst=FALSE
		}
	}
	
	if(postst && !kmst) {
		# boost provided startvalues using Viterbi algorithm
		if(class(dmm)[1]=="dmm") {
			logl=loglike(dat,dmm)$logl
			pp=dmm$pars
			postpars=poststart(dat,dmm)[1:dmm$npars]
			dmm$pars[which(dmm$fixed[1:dmm$npars]==1)]=postpars[which(dmm$fixed[1:dmm$npars]==1)]
			loglpost=loglike(dat,dmm)$logl
			if(loglpost<logl) dmm$pars=pp
			if(printlevel>9) cat("logl: ", logl, " post st logl: ", loglpost)
		}
	}
	
	xgmod=checkSetRecode(dat,dmm,tdcov,printlevel)
	
	if(tdcov) {
		if(der) {
			der=0
			warning("der set to FALSE for model with covariates.\n")
		}
		if(ses) {
			ses=FALSE
			warning("Can't compute standard errors with covariates model.\n Option ses set to FALSE\n")
		}
	}
	
	initLogl=loglike(dat=dat,dmm=xgmod,tdcov=tdcov,printlevel=0,set=FALSE)$logl
	if(is.nan(initLogl) || is.infinite(initLogl)) stop("Initial model inadmissible, better starting values needed?!\n")
	if(printlevel>0) cat("Initial loglikelihood: ", initLogl, " \n")
	
	timeUsed=0
	
	if(method=="npsol") {
		if(!is.loaded("npsolc")) {
			method="optim"
			warning("Optimization method changed to optim because npsol is not available on this computer.")
		}
	}
	
	# degrees of freedom in the data
	nobs1=0
	if(xgmod$ng==1) nobs1=sum(replicates(dat)*ntimes(dat))
	else for(i in 1:xgmod$ng) nobs1 = nobs1 + sum(replicates(dat[[i]])*ntimes(dat[[i]]))
	
## parameter reduction, selection etc	

	# model has no td pars or it has td pars and they should be fitted
	if(!xgmod$td || (xgmod$td==TRUE && tdcov==1)) {
		fitpars=xgmod$pars
		A=xgmod$A
		bl=xgmod$bl
		bllin=xgmod$bllin
		bu=xgmod$bu
		bulin=xgmod$bulin
		fixed=xgmod$fixed
	} else {
		# model has td pars but these should not be fitted
		# split pars in normal pars and tdpars
		fitpars=xgmod$pars[1:xgmod$npars]
		# split A, bl/u, and bl/ulin accordingly
		A=xgmod$A[,1:xgmod$npars,drop=FALSE]
		bl=xgmod$bl[1:xgmod$npars]
		bu=xgmod$bu[1:xgmod$npars]
		fixed=xgmod$fixed[1:xgmod$npars]
		
		# remove zero rows from A that result from removing tdpars, also remove rows with only one entry
		idx2=numeric(0)
		if(nrow(A)>0) {
			for(i in 1:nrow(A)) {
				if(sum(as.logical(A[i,]))==0 || sum(as.logical(A[i,]))==1) {
					idx2=c(idx2,i)
				}
			}
		}
		if(length(idx2)>0) {
			A=A[-idx2,,drop=FALSE]
			bllin=xgmod$bllin[-idx2]
			bulin=xgmod$bulin[-idx2]
		} else {
			bllin=xgmod$bllin
			bulin=xgmod$bulin
		}
	}
	
	# check for zero rows in A which may result from fixing parameters at zero, and remove them
	# set zeroes in A when pars are fixed, and change bllin and bulin accordingly
	# this is also done in function (mg-mix-)dmm, but users may add/change etc the linear constraint matrix 
	# after defining the model
	
	if(nrow(A)>0) {
		Bcomplete=matrix(as.logical(A),nrow=nrow(A))
		xcomplete=apply(Bcomplete,1,sum)
		A[,which(fixed==0)]=0
		B=matrix(as.logical(A),nrow=nrow(A))
		x=apply(B,1,sum)
		for(i in 1:nrow(A)) {
			if(x[i]!=xcomplete[i]) {
				for(j in 1:length(fitpars)) {
					if(fixed[j]==0 && Bcomplete[i,j]!=0 && fitpars[j]!=0) {
						bllin[i]=bllin[i]-bl[j]
						bulin[i]=bulin[i]-bu[j]
					}
				}
			}
		}		
		A=A[x>1,,drop=FALSE]
		bllin=bllin[x>1]
		bulin=bulin[x>1]
	}
	
	# only non-fixed pars are passed the optimization routine
	pars=fitpars[which(fixed==1)]
	
##  call npmain to optimize the model (non-linear constraints not implemented yet)
	if(method=="npsol") {
		
		bu=bu[which(fixed==1)]
		bl=bl[which(fixed==1)]
		
		if(printlevel>29) {
			print(A)
			print(bl)
			print(bu)
			print(bllin)
			print(bulin)
		}

		# npsol specific inputs
		fixedvals=fitpars
		A=A[,which(fixed==1),drop=FALSE]
		if(nrow(A)>0) {
			bl=c(bl,bllin)
			bu=c(bu,bulin)
		} # else nothing to do, bl and bu are already defined
		
		## npsol optimization options
#  		itlim=paste(c("Iteration Limit = ", iterlim , "   "), collapse=" ")
#  		npopts=c(itlim,
#  			"Print level = 20 ",
#  			"Minor print level = 10",
#  			"Derivative level = 0 ",
#  			"Hessian = No ",
#  			"Verify level = 0 ")
		
# 			"Optimality tolerance = 1.0e-12 ")
		# verify options: 0 is good for standard markovian and latent class models
		# verify level should be 3 for new types of models, ie with new distributions etc.
		
#  		optset=lapply(npopts,FUN=optfor <- function(opt) {.Fortran("npoptn",as.character(opt),PACKAGE="depmix") } )
		# optfile can be used to read in an optionsfile called npoptn from the working directory
		optfile=1 # using this will override above specified default options
		derivatives=der # are derivatives to be used or not
		cj=1 # non linear constraints are not supported yet, hence cj=1
		maxnpcalls=1 # maxnpcalls is not used at the moment
		nctotl=length(pars)+nrow(A)+0 # latter zero is for the number of non-linear constraints
		npars=length(pars) # nr of pars to be optimized
		
		timeUsed <-  system.time(z <- .C("npsolc",
				as.integer(npars),				# nr of pars to be optimized
				as.integer(nrow(A)), 			# nr of lin constraints
				as.integer(0), 					# nr of non-linear constraints
				as.integer(nrow(A)),			# lead dimension of A
				as.integer(1),					# lead dimension of Jacobian of non-linear constr
				as.integer(npars),				# lead dimension of matrix R
				as.double(A),					# linear constraint matrix
				as.double(bl),					# lower bounds on pars and constraints
				as.double(bu),					# upper bounds on pars and constraints
				inform=integer(1),				# inform
				iter=integer(1),				# nr of iterations
				istate=as.integer(rep(0,nctotl)),	# return value indicating whether constraints are satisfied
				as.double(1),					# par for non-linear constraints, not implemented
				as.double(0),					# not accessed, non-linear constraint var
				as.double(rep(0,nctotl)),		# clamda, bounds and linear constraints
				objf=double(1),					# return value: final logl
				gradu=double(npars),			# return value: gradients at final iterate
				R=as.double(rep(0,npars*npars)),# return value: augmented hessian
				pars=as.double(pars),			# initial values and return values with final values of pars
				totMem=integer(1),				# memory usage
				as.integer(maxnpcalls),			# may be used to restart iterations after perturbing parameters
				as.integer(optfile),			# logical indicating whether an options file should be read
				as.integer(printlevel),			# printlevel
				as.integer(derivatives),		# logical indicating whether analytical gradients are available or not
				as.integer(tdcov),				# logical indicating whether time dependent covariates pars are fitted or not
				as.integer(fixed),				# logical indicating which pars are fixed zeroes
				as.double(fixedvals),			# values of fixed parameters
				as.integer(length(fixed)),		# self evident
				PACKAGE="depmix"))
		z$timeUsed=timeUsed[3]
		if(printlevel>19) print(z)
		fitpars=xgmod$pars
		fitpars[which(fixed==1)]=z$pars
		z$pars=fitpars
		z
	}
	
 	##  call npmain to optimize the model (non-linear constraints not implemented yet)
 	if(method=="donlp") {
 		
 		requireNamespace("Rdonlp2")
 		
 		bu=bu[which(fixed==1)]
 		bl=bl[which(fixed==1)]
 		
 		if(printlevel>29) {
 			print(A)
 			print(bl)
 			print(bu)
 			print(bllin)
 			print(bulin)
 		}
 		
 		# donlp specific inputs
 		optpars=fitpars[which(fixed==1)]
 		A=A[,which(fixed==1),drop=FALSE]			
 		
 		# define loglike function
 		logl <- function(pars) {
 			xgmod$pars[which(fixed==1)]=pars
 			-loglike(dat=dat,dmm=xgmod,printlevel=0,set=FALSE,tdcov=tdcov)$logl
 		}
 		
 		if(der) {
 			grad <- function(pars) {
 				xgmod$pars[which(fixed==1)]=pars
 				gr <- -loglike(dat=dat,dmm=xgmod,printlevel=0,set=FALSE,tdcov=tdcov,grad=TRUE)$gr[which(fixed==1)]
 				return(gr)
			}
			
 			attr(logl, "gr") <- grad
 		}
 		
 		timeUsed <- system.time(
 			res <- Rdonlp2::donlp2(optpars, logl,
 				par.upper=bu,
 				par.lower=bl,
 				A = A,
 				lin.upper=bulin,
 				lin.lower=bllin,
 				nlin = list(),
 				control=Rdonlp2::donlp2Control(),
 				env=.GlobalEnv, name="Rdonlp2")
 		)
 		
 		z=list()
 		z$objf=-res$fx
 		z$iter=res$step.nr
 		z$inform=res$message
 		z$timeUsed=timeUsed[3]
 		fitpars[which(fixed==1)]=res$par
 		z$pars=fitpars
 		z$npars=xgmod$npars
 		z$totMem=NULL
 		z
 	}
	
	##  call optim or nlm to optimize the model
	if(method=="optim" || method=="nlm") {
# 		vfactor=15
		# reparametrize for optim and nlm
		noptpars=length(fitpars)
		# un-constrained pars
		uncon=rep(0,noptpars)
		constr=rep(0,noptpars)
		for(i in 1:noptpars) {
			if(fixed[i]==1) {
				if(sum(as.logical(xgmod$A[,i]))==0) uncon[i]=1
				else constr[i]=1
			}
		}
		# remove cols corresponding with non-constrained and fixed pars
		A=A[,which(constr==1),drop=FALSE]
		# remove rows correponding with inequality constraints
		idx=which(bllin<bulin) 
		ci=matrix(bulin,ncol=1) 
		
		ineq=FALSE
		if(length(idx)>0) {
			A=A[-idx,,drop=FALSE]
			ci=ci[-idx]
			Aineq=xgmod$A[idx,which(constr==1),drop=FALSE]
			ciIneqLow=matrix(xgmod$bllin[idx],ncol=1)
			ciIneqUp=matrix(xgmod$bulin[idx],ncol=1)
			ineq=TRUE
		}
		
		if(nrow(A)>0) {
			b=Null(t(A))
			ginvb=ginv(b)
			ginvA=ginv(A)
			newpars=ginvb%*%(fitpars[which(constr==1)]-ginvA%*%ci)
			nnewpars=length(newpars)
		} else {
			nnewpars=0
		}
		
		ucpars=xgmod$pars[which(uncon==1)]
		nucpars=length(ucpars)
		if(nnewpars>0) {
			pars=c(newpars,ucpars)
		} else {
			pars=ucpars
		}
		
		violationPenalty=vfactor
		
		viol <- function(xgpars) {
			violcon=0
			if(nrow(A)>0) violcon=violationPenalty*sum(abs(A%*%xgpars[which(constr==1)]-ci))
			violineq=0
			if(ineq) {
				violineqlow=violationPenalty*((Aineq%*%xgpars[which(constr==1)]<ciIneqLow))*abs(Aineq%*%xgpars[which(constr==1)]-ciIneqLow)
				violinequp=violationPenalty*((Aineq%*%xgpars[which(constr==1)]>ciIneqUp))*abs(Aineq%*%xgpars[which(constr==1)]-ciIneqUp)
				violineq=violineqlow+violinequp
			}
			violup=violationPenalty*sum(as.numeric(xgpars>xgmod$bu)*sqrt((xgpars-xgmod$bu)^2))
			viollow=violationPenalty*sum(as.numeric(xgpars<xgmod$bl)*sqrt((xgpars-xgmod$bl)^2))
			violation=violcon+violup+viollow+sum(violineq)
			return(violation)
		}
		
# 		assign("fcalls",0,envir=.GlobalEnv)
#   	assign("totalviolation",0,envir=.GlobalEnv)		
 		
		totalviolation = 0
		
		# optim optimization
		if(method=="optim") {
			f <- function(pars) {
				if(nnewpars>0) xgmod$pars[which(constr==1)]=ginvA%*%ci + b%*%pars[1:nnewpars]
				if(nucpars>0) xgmod$pars[which(uncon==1)]=pars[(nnewpars+1):(nnewpars+nucpars)]
				ll=loglike(dat=dat,dmm=xgmod,printlevel=0,set=FALSE,tdcov=tdcov)$logl
				xgpars=xgmod$pars
				violation=viol(xgpars)
#  				assign("totalviolation",violation,envir=.GlobalEnv)
# 				totalviolation=violation
				totalviolation=viol(xgpars)
# 				assign("fcalls",fcalls+1,envir=.GlobalEnv)
				if(is.nan(ll) || is.infinite(ll)) ll=initLogl*1.01
				ll=ll-totalviolation
				pr=0
# 				if(printlevel>9 && (fcalls%%10)==0) {cat("call " , fcalls, " logl ", ll, "\n"); pr=1}
# 				if(printlevel>14 && pr==0) cat("call " , fcalls, " logl ", ll, "\n")
				ll
			}
			
			gr <- function(pars) {
				if(nnewpars>0) xgmod$pars[which(constr==1)]=ginvA%*%ci + b%*%pars[1:nnewpars]
				if(nucpars>0) xgmod$pars[which(uncon==1)]=pars[(nnewpars+1):(nnewpars+nucpars)]
				gr=loglike(dat=dat,dmm=xgmod,grad=TRUE,printlevel=0,set=FALSE)$gr
				if(nnewpars>0) grad=c(gr[which(constr==1)]%*%b,gr[which(uncon==1)])
				else grad=gr[which(uncon==1)]
# 				if(printlevel>0 && (fcalls%%10)==0) cat("call " , fcalls, " grad ", grad, "\n")
# 				if(printlevel>9) cat("call " , fcalls, " grad ", grad, "\n")
				grad
			}
			
			# optim control parameters
			nopt=1
			fnscale=-1.0
			trace=1
			factr=1e7
			maxit=5000
			REPORT=10
			if(printlevel>1 && printlevel<5) REPORT=5
			if(printlevel>4) REPORT=1
			if(printlevel==0) trace=0
			control=list(fnscale=fnscale,trace=trace,factr=factr,maxit=maxit,REPORT=REPORT)
			oldObjf=initLogl
			
			for(i in 1:nopt) { 
				timeUsed=timeUsed+system.time(z<-optim(pars,f,gr=NULL,control=control,method="BFGS",hessian=FALSE))[3]
				if(nnewpars>0) xgmod$pars[which(constr==1)]=ginvA%*%ci + b%*%z$par[1:nnewpars]
				if(nucpars>0) xgmod$pars[which(uncon==1)]=z$par[(nnewpars+1):(nnewpars+nucpars)]
				
				if(totalviolation > 1e-9) {
					cat("Constraint violation not zero, minimum probably not optimal, provide better starting values.\n")
					cat("Violation is: ", totalviolation, "\n")
				}
				
				if(z$convergence==0 && totalviolation < 1e-9) break
				violationPenalty=3*violationPenalty
				oldObjf=z$value
				pars=z$par
				
			}
			names(z)=c("par","objf","iter","inform","message")
		} # end optim
		
		# nlm optimization
		if(method=="nlm") {
# 			assign("loglvalue",0,envir=.GlobalEnv)
			
			f <- function(pars) {
				# expand pars into the model
				if(nnewpars>0) xgmod$pars[which(constr==1)]=ginvA%*%ci + b%*%pars[1:nnewpars]
				if(nucpars>0) xgmod$pars[which(uncon==1)]=pars[(nnewpars+1):(nnewpars+nucpars)]
				# get the actual likelihood and gradients if requested/feasible
				if(der) {
					llgr=loglike(dat=dat,dmm=xgmod,printlevel=0,grad=TRUE,set=FALSE)
					if(nnewpars>0) grad=c(-llgr$gr[which(constr==1)]%*%b,-llgr$gr[which(uncon==1)])
					else grad=-llgr$gr[which(uncon==1)]
				} else llgr=loglike(dat=dat,dmm=xgmod,tdcov=tdcov,printlevel=0,grad=FALSE,set=FALSE)
				ll=llgr$logl
				violation=viol(xgmod$pars)
				totalviolation=viol(xgmod$pars)
# 				totalviolation=violoation
# 				assign("totalviolation",violation,envir=.GlobalEnv)
				ll=ll-totalviolation
				pr=0
# 				assign("fcalls",fcalls+1,envir=.GlobalEnv)
				if(is.nan(ll) || is.infinite(ll)) ll=initLogl*1.01
# 				if(printlevel>9 && (fcalls%%10)==0) {cat("call " , fcalls, " logl ", ll, "\n"); pr=1}
# 				if(printlevel>14 && pr==0) cat("call " , fcalls, " logl ", ll, "\n")
				if(der) attr(ll,'gradient')=grad
				-ll
			}
			
			prl=1
			if(printlevel==0) prl=0
			if(printlevel>4) prl=2
			ndig=12
			nopt=1
			adddig=3
			
			for(i in 1:nopt) {
				gtol=10^-(ndig-6)
				stol=10^-(ndig-6)
				timeUsed=timeUsed+system.time(z<-nlm(f,p=pars,hessian=FALSE,ndigit=ndig,gradtol=gtol,steptol=stol,fscale=-initLogl*0.95,print.level=prl,check.analyticals = FALSE))[3]
				if(nnewpars>0) xgmod$pars[which(constr==1)]=ginvA%*%ci + b%*%z$estimate[1:nnewpars]
				if(nucpars>0) xgmod$pars[which(uncon==1)]=z$estimate[(nnewpars+1):(nnewpars+nucpars)]
				if(z$code<3 && totalviolation<1e-9) break
				if(z$code==3) { ndig=ndig+adddig; vfactor=3*vfactor}
				if(totalviolation > 1e-9) vfactor=3*vfactor
				if(z$code>3) break
			}
			
			if(totalviolation > 1e-9) {
				cat("Constraints not satisfies provide better starting/fixed values or increase vfactor.\n")
				cat("Violation is: ", totalviolation, "\n")
			}
			
			names(z)=c("objf","par","grad","inform","iter")
# 			z$objf=ifelse(loglvalue==0,-z$objf,loglvalue)
		} # end nlm
		
		# put fitted pars into xgmod, gather output into z
		if(nnewpars>0) xgmod$pars[which(constr==1)]=ginvA%*%ci + b%*%z$par[1:nnewpars]
		if(nucpars>0) xgmod$pars[which(uncon==1)]=z$par[(nnewpars+1):(nnewpars+nucpars)]
			
		z$timeUsed=timeUsed
		z$pars=xgmod$pars
		z$npars=xgmod$npars
		z$totMem=NULL
	} # end optim/nlm optimization
	
	if(printlevel>0) cat("Final loglikelihood: ", z$objf,"\n")
	
# put fitted pars into the original model
	mod$pars=z$pars
	if(class(mod)[1]=="mixdmm") {
		bp=mod$nrcomp
		for(i in 1:mod$nrcomp) {
			mod$mod[[i]]$pars = c(z$pars[i],z$pars[(bp+1):(bp+mod$mod[[i]]$npars-1)])
			if(mod$mod[[i]]$td) {
				if(tdcov) {
					mod$mod[[i]]$pars = c(mod$mod[[i]]$pars,z$pars[i+mod$npars],z$pars[(bp+1+mod$npars):(bp+mod$mod[[i]]$npars-1+mod$npars)])
					mod$mod[[i]]$tdfit=1
				} else {
					mod$mod[[i]]$tdfit=0
				}
			}
			bp = bp + mod$mod[[i]]$npars-1
		}
	}
	if(class(mod)[1]=="mgd") {
		pergrouppars = mod$npars/mod$ng
		bg=0
		for(g in 1:mod$ng) {
			mixpars=z$pars[(bg+1):(bg+pergrouppars)]
			if(mod$td) mixtdpars=z$pars[(bg+1+mod$npars):(bg+pergrouppars+mod$npars)]
			mod$mixmod[[g]]$pars=mixpars
			bp=mod$nrcomp
			if(bp>1) {
				for(i in 1:mod$nrcomp) {
					mod$mixmod[[g]]$mod[[i]]$pars = c(mixpars[i],mixpars[(bp+1):(bp+mod$mixmod[[g]]$mod[[i]]$npars-1)])
					if(mod$mixmod[[g]]$mod[[i]]$td) {
						if(tdcov) {
							mod$mixmod[[g]]$mod[[i]]$pars = c(mod$mixmod[[g]]$mod[[i]]$pars,mixtdpars[i],mixtdpars[(bp+1):(bp+mod$mixmod[[g]]$mod[[i]]$npars-1)])
							mod$mixmod[[g]]$mod[[i]]$tdfit=1
						} else {
							mod$mixmod[[g]]$mod[[i]]$tdfit=0
						}
					}
					bp = bp + mod$mixmod[[g]]$mod[[i]]$npars-1
				}
			} else {
				mod$mixmod[[g]]$pars = c(mixpars[1],mixpars[(bp+1):(bp+mod$mixmod[[g]]$npars-1)])
				if(mod$mixmod[[g]]$td) {
					if(tdcov) {
						mod$mixmod[[g]]$pars = c(mod$mixmod[[g]]$pars,mixtdpars[1],mixtdpars[(bp+1):(bp+mod$mixmod[[g]]$npars-1)])
						mod$mixmod[[g]]$tdfit=1
					} else {
						mod$mixmod[[g]]$tdfit=0
					}
				}
			}
			bg=bg+pergrouppars
		}
	}
	if(mod$td && tdcov==0) mod$tdfit=0
	if(mod$td && tdcov==1) mod$tdfit=1
		
## compute posteriors
	if(poster) {
		if(printlevel>0) cat("Computing posteriors \n")
		x<-system.time(pp<-posterior(dat=dat,dmm=mod,tdcov,printlevel))
		z$timeUsed=z$timeUsed+x[3]
	}
	if(ses) {
		val=computeSes(dat,mod)
		se=val$se
		hs=val$hs
	}
	
## fit measures
	if(!tdcov && mod$td) {
		df=mod$freeparsnotd
	} else df=mod$freepars
	aic = -2*z$objf+2*df
	bic = -2*z$objf+log(nobs1)*df
	
	if(printlevel>0) {
		cat("This took ", sum(z$iter), "iterations, ", z$timeUsed," seconds \n")
		if(method=="npsol") cat(" and ", z$totMem, " bytes of memory","\n")
		if(method=="optim" && z$inform>0) cat(z$message,"\n")
	}
	
## return values
	fit=list(date=date(),timeUsed=z$timeUsed,totMem=z$totMem,iter=z$iter,inform=z$inform,
	loglike=z$objf,aic=aic,bic=bic,df=df,nobs1=nobs1,mod=mod)
	if(poster) fit$post=pp
	if(ses) {fit$se=se; fit$hs=hs}
	fit$method=method
	class(fit)="fit"
	return(fit)
}

loglike <- function(dat,dmm,tdcov=0,grad=FALSE,hess=FALSE,set=TRUE,grInd=0,sca=1.0,printlevel=1) {
	
	if(class(dmm)[1]=="fit") dmm <- dmm$mod
	
	if(set) {
		xgmod=checkSetRecode(dat,dmm,tdcov,printlevel)
	} else {
		if(class(dmm)[1]!="mgd") xgmod <- mgdmm(dmm=dmm)
		else xgmod <- dmm
	}
	
	if(any(xgmod$itemtypes< 0 )) {
		grad=0
		warning("gradients not available for items other than multinomial and gaussian.")
	}
	
	if(hess) {
		warning("Hessian not implemented (yet).\n")
		hess=0
	}
	
	if(tdcov) {
		pars=xgmod$pars
		grad=0
		hess=0
		gr=1
		hs=1
		grset=1
		hsset=1
		grInd=0
		nrGradInd=1
	} else {
		pars=xgmod$pars[1:xgmod$npars]
		if(grad) gr=length(pars)
		else gr=1
		if(hess) hs=length(pars)*length(pars)
		else hs=1
		if(xgmod$ng==1 && xgmod$nrcomp==1) {
			grInd=grInd
			nrGradInd=length(ntimes(dat))*xgmod$npars
		} else {
			grInd=0
			nrGradInd=1
		}
	}
	
	z <- .C("loglikelihood",
		as.double(pars),
		as.integer(tdcov),
		objf=as.double(0.0),
		grad=as.integer(grad),
		hess=as.integer(hess),
		gr=double(gr),
		hs=double(hs),
		grset=integer(gr),
		hsset=integer(hs),
		as.integer(grInd),
		gradIndividual=double(nrGradInd),
		sc=as.double(sca),
		info=integer(1),
		PACKAGE="depmix")
	
	value=list(logl=z$objf)
	if(grad) {
		value$gr=z$gr
		value$grset=z$grset
	}
	if(hess) {
		value$hs=matrix(z$hs,length(pars),length(pars))
		value$hsset=z$hsset
	}
	if(grInd) value$gradIndividual=z$gradIndividual
	value
}

# compute posterior probabilities
posterior <- function(dat,dmm,tdcov=0,printlevel=1) {
	
	if(class(dmm)[1]=="fit") dmm <- dmm$mod
	
	xgmod=checkSetRecode(dat,dmm,tdcov,printlevel)
	if(xgmod$ng==1) dat=list(dat)
	
	post <- list()
	post$states <- list() 
	post$comp <- list() 
	for(g in 1:xgmod$ng) {
		post$comp[[g]]=matrix(0,length(ntimes(dat[[g]])),xgmod$nrcomp+1)
		post$states[[g]] <- matrix(0,ncol=(2+sum(xgmod$nstates)),nrow=0)
	###the first entry in each series is the number of the component
		for(i in 1:(length(ntimes(dat[[g]])))) {
			pp=ntimes(dat[[g]])[i]*sum(xgmod$nstates)
			z <- .C("posteriors",
				states=double(ntimes(dat[[g]])[i]),
				postdelta=double(pp),
				postcomp=integer(1),
				as.double(xgmod$pars),
				as.integer(tdcov),
				as.integer(g),
				as.integer(i),
				PACKAGE="depmix")
			post$comp[[g]][i,xgmod$nrcomp+1]=z$postcomp
			pst=matrix(c(rep(z$postcomp,ntimes(dat[[g]])[i]),z$states,matrix(z$postdelta,ncol=sum(xgmod$nstates),byrow=TRUE)),ncol=(2+sum(xgmod$nstates)))
			postprob=pst[ntimes(dat[[g]])[i],3:(2+sum(xgmod$nstates))]
			bg=1; en=xgmod$nstates[1];
			post$comp[[g]][i,1] = sum(postprob[bg:en])
			if(xgmod$nrcomp>1) {
				for(cp in 2:xgmod$nrcomp) {	
					bg=en+1; en=bg+xgmod$nstates[cp]-1
					post$comp[[g]][i,cp] = sum(postprob[bg:en])
				}
			}
			post$states[[g]]=rbind(post$states[[g]],pst)
		}
	}
	post
}

computeSes <- function(dat,dmm) {
	
	if(class(dmm)[1]=="fit") dmm <- dmm$mod

	cat("Computing standard errors \n")
	mod=dmm
	eps=1e-8
	lgh=loglike(dat=dat,dmm=mod,grad=TRUE,tdcov=0,printlevel=0)
	ll=lgh$logl
	gr=lgh$gr
	hsall=matrix(0,mod$npars,mod$npars)
	for(j in 1:mod$npars) {
		if(mod$fixed[j]!=0) {
			del=eps*max(1,mod$pars[j])
			por=mod$pars[j]
			mod$pars[j] = por+del
			alli=loglike(dat=dat,dmm=mod,grad=TRUE,hess=FALSE,set=FALSE,tdcov=0,printlevel=0)
			grd=alli$gr
			for(i in 1:mod$npars) hsall[i,j] = (grd[i]-gr[i])/del
			mod$pars[j]=por
		}
	}
	se=rep(0,mod$npars)
	hs=-hsall[which(mod$fixed[1:mod$npars]==1),which(mod$fixed[1:mod$npars]==1)]	
	if(nrow(mod$A)>0) {
		A=matrix(mod$A[,which(mod$fixed[1:mod$npars]==1)],nrow=nrow(mod$A))
		idx=which(mod$bllin<mod$bulin)
		if(length(idx)>0) A=A[-idx,,drop=FALSE]
		idx=numeric(0)
		for(i in 1:nrow(A)) {
			if(sum(as.logical(A[i,]))==0) {
				idx=c(idx,i)
			}
		}
		if(length(idx)>0) A=A[-idx,,drop=FALSE]
		d=hs+t(A)%*%A
		di=try(solve(d),silent=TRUE)
		if(class(di)=="try-error") {
			warning("Hessian singular, ses could not be computed.") 
			val=list(se=0,hs=0)
		} else {
			ada=A%*%di%*%t(A)
			adai=try(solve(ada),silent=TRUE)
			if(class(adai)=="try-error") {
				warning("Near-singular hessian, ses may be bad.\n")
				diag(ada)=diag(ada)*1.000001
				adai=try(solve(ada))
				if(class(adai)=="try-error") {
					warning("Corrected hessian also singular, ses computed without contraints.\n")
					se[which(mod$fixed[1:mod$npars]==1)]=sqrt(diag(di))
				} else {
					ch=di-di%*%t(A)%*%adai%*%A%*%di
					se[which(mod$fixed[1:mod$npars]==1)]=sqrt(diag(ch))
				}
			} else {
				ch=di-di%*%t(A)%*%adai%*%A%*%di
				se[which(mod$fixed[1:mod$npars]==1)]=sqrt(diag(ch))
			} 
		} 
	} else {
		se[which(mod$fixed[1:mod$npars]==1)]=sqrt(diag(solve(hs)))
	}
	val=list(se=se,hs=hs)
	val
	
}

# compute bootstrap based standard errors + gof
bootstrap <- function(object, dat, samples=100, pvalonly=0, ...) {
	if(!class(object)=="fit") stop("Bootstrapping only possible on 'fit' objects.")
	mod=object$mod
	# compute only bootrapped pval, not ses
	gof=matrix(0,samples,3)
	bootpars=matrix(0,samples,mod$npars)
	llcrit=0
	if(pvalonly==1) {
		llcrit=object$loglike-0.00001*object$loglike
		.C("setCrit",
			as.double(llcrit),
			PACKAGE="depmix")
	}
	for(i in 1:samples) {
		gen=generate(dmm=mod, ntimes=ntimes(dat))
		ll=loglike(dat=gen,dmm=mod,printlevel=0)$logl
		fit=fitdmm(dat=gen,dmm=mod,printlevel=0,ses=0,postst=0,poster=FALSE)
		bootpars[i,]=fit$mod$pars
		ll=fit$loglike
		gof[i,]=c(ll,fit$aic,fit$bic)
		cat(i, " ", ll, "\n")
	}
	better=matrix(0,samples,1)
	better=ifelse(gof[,1]>object$logl,1,0)
	if(pvalonly==0) {
		object$bootpars=bootpars
		object$bse=apply(bootpars,2,sd)
		object$gof=gof
	}
	else object$gof=gof[,1]
	object$better=better
	object$pbetter=sum(better)/samples
	object
}

summary.fit <- function(object,precision=3, fd=1, ...) {
	cat("Model: ", object$mod$modname," fitted at ", object$date, " \n")
	cat("Optimization information, method is ", object$method, "\n")
	cat(" Iterations: ", object$iter)
	cat("\n Inform: ", object$inform, " (look up the respective manuals for more information.)\n")
	cat("\n Loglikelihood of fitted model: ", round(object$loglike,precision))
	cat("\n AIC: ",  round(object$aic,precision))
	cat("\n BIC: ",  round(object$bic,precision))
	cat("\n Number of observations (used in BIC): ", object$nobs)
	cat("\n Fitted model \n")
	if(!is.null(object$bse) && !is.null(object$se)) {
		if(fd==0) summary(object$mod,precision=precision, se=object$bse)
		else summary(object$mod,precision=precision, se=object$se)
	} else {
		if(!is.null(object$se)) summary(object$mod,precision=precision, se=object$se)
		if(!is.null(object$bse)) summary(object$mod,precision=precision, se=object$bse)
	}
	if(is.null(object$se) && is.null(object$bse)) summary(object$mod, precision=precision)
}

oneliner=function(object,precision=3) {
	x=c(object$mod$modname,round(object$loglike,precision),round(object$aic,precision),round(object$bic,precision),
		object$mod$npars,object$mod$freepars,object$date)
	return(x)
}

