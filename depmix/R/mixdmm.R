
########################
#                      #
#   MIXTURE OF DMM'S   #
#                      #
########################

mixdmm <- function(dmm, modname=NULL, mixprop=NULL, conrows=NULL) {
	nrcomp=length(dmm)
	for(i in 1:nrcomp) if(!(class(dmm[[i]])[1]=="dmm")) stop("Argument dmm must be a list of models with class dmm.")
	if(!(is.null(conrows))) {
		if((length(conrows)%%nrcomp)!=0) stop("Length(conrows) is not a multiple of the number of components")}
	if(is.null(mixprop)) mixprop=rep(1/nrcomp,nrcomp)
	else if(length(mixprop)!=nrcomp) stop("mixprop has incorrect length")
	mixprop <- mixprop/sum(mixprop)
	if(is.null(modname)) modname=paste(nrcomp,"-component mixture model")
	nstates=numeric(0)
	mod=dmm
	nstates = mod[[1]]$nstates
	if(nrcomp>1) for(i in 2:nrcomp) nstates=c(nstates,mod[[i]]$nstates)
	itemtypes=dmm[[1]]$itemtypes
# single component model, not much to do
	if(nrcomp==1) {
		mixmod=mod[[1]]
		mixmod$nrcomp=1
		if(mod[[1]]$td) td=TRUE
		else td=FALSE
		class(mixmod)="dmm"
	} else { 
# nrcomp>1
# paramaters
	# bounds on mixprop parameters
		bl=rep(0,nrcomp)
		bu=rep(1,nrcomp)
	# the default is that mixprop pars are free
		fixed=rep(1,nrcomp)
	# counting parameters
		npars = 0
	# join parameters, bounds and fixed values
		pars=mixprop
		td=FALSE
		for(i in 1:nrcomp) {
			npars = npars+mod[[i]]$npars
			if(mod[[i]]$td) td=TRUE
			pars = c(pars,mod[[i]]$pars[2:mod[[i]]$npars])
			bl=c(bl,mod[[i]]$bl[2:mod[[i]]$npars])
			bu=c(bu,mod[[i]]$bu[2:mod[[i]]$npars])
			fixed=c(fixed,mod[[i]]$fixed[2:mod[[i]]$npars])
		}
# sum constraint for mixprop
		sc=matrix(rep(1,nrcomp),nrow=1)
		bllin=1
		bulin=1
	# user supplied general linear constraints on mixprop
		amp=NULL
		if(!(is.null(conrows))) {
			conrows=matrix(conrows,ncol=npars,byrow=TRUE)
			amp=conrows
			bllin=c(bllin,rep(0,nrow(amp)))
			bulin=c(bulin,rep(0,nrow(amp)))
		}
		if(is.null(amp)) amp=sc 
		else amp=rbind(sc,amp)
# join constraint matrices
		nr=nrow(amp)
	# counting constraints and setting their bounds
		linmat=list()
		for(i in 1:nrcomp) {
			if(nrow(mod[[i]]$A)>0) {
				nr = nr + nrow(mod[[i]]$A)
				# this is assuming there are no constraints between components
				bllin=c(bllin,mod[[i]]$bllin)
				bulin=c(bulin,mod[[i]]$bulin)
				linmat[[i]]=mod[[i]]$A[,2:mod[[i]]$npars,drop=FALSE]
			} else {
				# we need a matrix with positive number of rows here, to work in function bdiag, it will 
				# automatically removed again later
				linmat[[i]]=matrix(0,1,mod[[i]]$npars-1)
				bllin=c(bllin,0)
				bulin=c(bulin,0)
				nr=nr+1
			}
		}
# make combined constraint matrix 
		A=matrix(0,ncol=npars,nrow=nr)
		A[1:nrow(amp),1:nrcomp] <- amp
		B=bdiag(linmat)
		A[2:nrow(A),(nrcomp+1):npars]=B
		
		freeparsnotd=sum(as.logical(fixed))
		Ared=A[,which(fixed==1),drop=FALSE]
		freeparsnotd=freeparsnotd-qr(Ared[which(bllin==bulin),,drop=FALSE])$rank
		
# covariate constraints and parameters
		nparstotal=npars
		if(td) {
			# add zeroes for all possible pars, it is possible that some 
			# components have no td covariates
			pars = c(pars,rep(0,npars))
			bl=c(bl,rep(0,npars))
			bu=c(bu,rep(0,npars))
			fixed=c(fixed,rep(0,npars))
			linmattd=list()
			for(i in 1:nrcomp) {
				# if the component has covariates, add the appropriate parameter values, 
				# bounds and fixed constraints 
				if(mod[[i]]$td) { ###										where does this 2 come from??? ****
 					pars[paridx(nstates,itemtypes,mat="all",comp=i)+npars]=mod[[i]]$pars[(mod[[i]]$npars+2):(2*mod[[i]]$npars)]
					bl[paridx(nstates,itemtypes,mat="all",comp=i)+npars]=mod[[i]]$bl[(mod[[i]]$npars+2):(2*mod[[i]]$npars)]
					bu[paridx(nstates,itemtypes,mat="all",comp=i)+npars]=mod[[i]]$bu[(mod[[i]]$npars+2):(2*mod[[i]]$npars)]
					fixed[paridx(nstates,itemtypes,mat="all",comp=i)+npars]=mod[[i]]$fixed[(mod[[i]]$npars+2):(2*mod[[i]]$npars)]
				}
				if(nrow(mod[[i]]$A)>0) {
					if(mod[[i]]$td) linmattd[[i]]=mod[[i]]$A[,(mod[[i]]$npars+2):(mod[[i]]$nparstotal)]
					else linmattd[[i]]=matrix(0,nrow(mod[[i]]$A),mod[[i]]$npars-1)
				} else {
					# we need a matrix with positive number of rows here, to work in function bdiag, it will 
					# automatically removed again later
					linmattd[[i]]=matrix(0,1,mod[[i]]$npars-1)
				}
			}
			zeroes=matrix(0,ncol=npars,nrow=nr)
			A=cbind(A,zeroes)
			C=bdiag(linmattd)
			A[2:nrow(A),(npars+nrcomp+1):(2*npars)]=C
			nparstotal=2*npars
		}
		for(i in 1:nparstotal) { 
			if(fixed[i]==0) { 
				bl[i]=pars[i] 
				bu[i]=pars[i] 
			}
		}
		for(i in 1:nparstotal) {if(fixed[i]==0) A[,i]=0}
		freepars=sum(as.logical(fixed))-qr(A[which(bllin==bulin),])$rank
		mixmod=list(nrcomp=nrcomp,nstates=nstates,itemtypes=itemtypes,mod=mod,modname=modname,
			npars=npars,nparstotal=nparstotal,freepars=freepars,
			freeparsnotd=freeparsnotd,pars=pars,fixed=fixed,
			A=A,bl=bl,bu=bu,bllin=bllin,bulin=bulin,td=td,tdfit=1,st=mod[[1]]$st)
		if(td) mixmod$td=TRUE
		class(mixmod) = "mixdmm"
	} # end nrcomp>1
	mixmod
}

summary.mixdmm <- function(object, specs=FALSE, precision=3, se=NULL, ...) {
	cat(" Model: ", object$modname,"\n")
	cat(" Nr of components: ", object$nrcomp,"\n")
	cat(" npars:            ", object$nparstotal, "\n")
	cat(" freepars:         ", ifelse(object$tdfit,object$freepars,object$freeparsnotd), "\n\n")
	cat(" Mixture proportions \n")
	mpv=object$pars[1:object$nrcomp]
	ses=FALSE
	if(!is.null(se)) {
		ses=TRUE
		mpse=se[1:object$nrcomp]
		mpt=rep(0,object$nrcomp)
		for(i in 1:object$nrcomp) {
			if(mpse[i]!=0) mpt[i]=mpv[i]/mpse[i]
			else mpt[i]=NA
		}
		se=se[(object$nrcomp+1):object$npars]
	}
	if(ses) {
		mp=matrix(c(mpv,mpse,mpt),3,object$nrcomp,byrow=TRUE)
		rownames(mp)=c("val","se","t") 
		oln=rep("Comp ",object$nrcomp)
		oln=paste(oln,1:object$nrcomp)
		colnames(mp)=oln
		print(mp)
	}
	else print(mpv)
	for(i in 1:object$nrcomp) {
		cat("\n Component: ", i, "\n")
		if(ses) {
			nrs=object$mod[[i]]$npars-1
			summary(object$mod[[i]],specs,precision, se=c(0,se[1:nrs]), ...)
			se=se[(nrs+1):length(se)]
		} else summary(object$mod[[i]],specs,precision, ...)
	}
}
