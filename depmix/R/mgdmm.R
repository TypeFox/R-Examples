
#########################
#                       #
#   MULTI GROUP MODEL   #
#                       #
#########################

mgdmm <- function(dmm,ng=1,modname=NULL,trans=FALSE,obser=FALSE,init=FALSE,conpat=NULL) {
# create a (mixture) model for each group, if not provided
	mixmod=list()
	# separate models provided for each group
	if(length(dmm)==ng) {
		if(class(dmm[[1]])[1]=="dmm") {
			for(i in 1:ng) {
				mixmod[[i]]=mixdmm(dmm=list(dmm[[i]]))
			}
		} else {
			mixmod <- dmm
		}
	} else {
	# each group gets the same model, ie with the same constraints, start values etc.
		if(class(dmm)[1]=="dmm") mixmods=mixdmm(dmm=list(dmm))
		else mixmods=dmm
		for(i in 1:ng) mixmod[[i]] <- mixmods
	}
	xmod=mixmod[[1]]
	nrcomp=xmod$nrcomp
# ng=1, not much to do
	if(ng==1) {
		mgd=mixmod[[1]]
		if(mixmod[[1]]$td) td=TRUE
		else td=FALSE
		mgd$ng=1
	} else { # ng>1
# count parameters etc.
		nstates=xmod$nstates
		itemtypes=xmod$itemtypes
		npars=xmod$npars*ng
# join parameters, bu, bl and fixed
		pars=numeric(0); bl=numeric(0); bu=numeric(0); fixed=numeric(0)
		linmat=list()
		bllin=numeric(0); bulin=numeric(0)
		td=FALSE
		for(i in 1:ng) {
			if(mixmod[[i]]$td) td=TRUE
			pars=c(pars,mixmod[[i]]$pars[1:mixmod[[i]]$npars])
			bl=c(bl,mixmod[[i]]$bl[1:mixmod[[i]]$npars])
			bu=c(bu,mixmod[[i]]$bu[1:mixmod[[i]]$npars])
			fixed=c(fixed,mixmod[[i]]$fixed[1:mixmod[[i]]$npars])
			linmat[[i]]=mixmod[[i]]$A[,1:mixmod[[i]]$npars,drop=FALSE]
			bllin=c(bllin,mixmod[[i]]$bllin)
			bulin=c(bulin,mixmod[[i]]$bulin)
		}
# put A together for multiple groups
 		A=bdiag(linmat)
		nr=nrow(A)
# between group contraints
	# trans parameters equal
		if(trans) {
			# count the number of transition pars and add this many rows *(ng-1) to A
			trcount=0
			for(i in 1:nrcomp) {
				trcount = length(paridx(nstates,itemtypes,comp=i,mat="tr"))
				trinv = matrix(0,trcount*(ng-1),npars)
				j=0
				for(g in 1:(ng-1)) {
					idx = paridx(nstates,itemtypes,comp=i,mat="tr")
					for (id in idx) {
						j = j+1
						trinv[j,id] = 1
						trinv[j,id+xmod$npars*g]=-1
					}
				}				
				A <- rbind(A,trinv)
				bltr = rep(0,trcount*(ng-1))
				butr = rep(0,trcount*(ng-1))
				bllin = c(bllin,bltr)
				bulin = c(bulin,butr)
			}
		} # end trans
		# obser parameters equal
		if(obser) {
			# count the number of observation pars and add this many rows *(ng-1) to A
			obcount=0
			for(i in 1:nrcomp) {
				obcount = length(paridx(nstates,itemtypes,comp=i,mat="ob"))
				obinv = matrix(0,obcount*(ng-1),npars)
				j=0
				for(g in 1:(ng-1)) {
					idx = paridx(nstates,itemtypes,comp=i,mat="ob")
					for(id in idx) {
						j = j+1
						obinv[j,id] = 1
						obinv[j,id+xmod$npars*g]=-1
					}
				}
				A <- rbind(A,obinv)
				blob = rep(0,obcount*(ng-1))
				buob = rep(0,obcount*(ng-1))
				bllin = c(bllin,blob)
				bulin = c(bulin,buob)
			}
		}
		# init parameters equal
		if(init) {
			# count the number of initial state pars and add this many rows *(ng-1) to A
			incount=0
			for(i in 1:nrcomp) {
				incount = length(paridx(nstates,itemtypes,comp=i,mat="in"))
				ininv = matrix(0,incount*(ng-1),npars)
				j=0
				for(g in 1:(ng-1)) {
					idx = paridx(nstates,itemtypes,comp=i,mat="in")
					for (id in idx) {
						j = j+1
						ininv[j,id] = 1
						ininv[j,id+xmod$npars*g]=-1
					}
				}
				A <- rbind(A,ininv)
				blin = rep(0,incount*(ng-1))
				buin = rep(0,incount*(ng-1))
				bllin = c(bllin,blin)
				bulin = c(bulin,buin)
			}
		}
		nrinv=nrow(A)-nr
		freeparsnotd=sum(as.logical(fixed))
		Ared=A[,which(fixed==1),drop=FALSE]
		freeparsnotd=freeparsnotd-qr(Ared[which(bllin==bulin),,drop=FALSE])$rank
# covariate constraints and stuff like that
		nparstotal=npars
		if(td) {
			# add zeroes for all possible pars, it is possible that some 
			# components or groups have no td covariates
			pars = c(pars,rep(0,npars))
			bl=c(bl,rep(0,npars))
			bu=c(bu,rep(0,npars))
			fixed=c(fixed,rep(0,npars))
			nrcp=rep(0,nrcomp)
			if(nrcomp>1) {
				for(i in 1:nrcomp) {
					nrcp[i]=mixmod[[1]]$mod[[i]]$npars
				}
			} else nrcp=mixmod[[1]]$npars
			ngpars=sum(nrcp)
			linmattd=list()
			for(i in 1:ng) {
				# if the component has covariates, add the appropriate parameter values, 
				# bounds and fixed constraints 
				if(mixmod[[i]]$td) {
					pars[(1+npars):(mixmod[[i]]$npars+npars+ngpars*(i-1))]=mixmod[[i]]$pars[(mixmod[[i]]$npars+1):(mixmod[[i]]$nparstotal)]
					bl[(1+npars):(mixmod[[i]]$npars+npars+ngpars*(i-1))]=mixmod[[i]]$bl[(mixmod[[i]]$npars+1):(mixmod[[i]]$nparstotal)]
					bu[(1+npars):(mixmod[[i]]$npars+npars+ngpars*(i-1))]=mixmod[[i]]$bu[(mixmod[[i]]$npars+1):(mixmod[[i]]$nparstotal)]
					fixed[(1+npars):(mixmod[[i]]$npars+npars+ngpars*(i-1))]=mixmod[[i]]$fixed[(mixmod[[i]]$npars+1):(mixmod[[i]]$nparstotal)]
				}
 				if(nrow(mixmod[[i]]$A)>0) {
 					if(mixmod[[i]]$td) linmattd[[i]]=mixmod[[i]]$A[,(mixmod[[i]]$npars+1):(mixmod[[i]]$nparstotal)]
 					else linmattd[[i]]=matrix(0,nrow(mixmod[[i]]$A),mixmod[[i]]$npars)
 				} else {
 					linmattd[[i]]=matrix(0,0,mixmod[[i]]$npars)
 				}
			}
			nparstotal=npars*2
			Atd=bdiag(linmattd)
			Atd=rbind(Atd,matrix(0,nrinv,npars))
 			A=cbind(A,Atd)
		} # end td pars
		if(!td && !is.null(conpat)) {
			fixed=pa2conr(conpat)$fix
			conp=pa2conr(conpat)$conr #pa2conr in depmix-internal
			if(nrow(conp)>0) {
				A=rbind(A,conp)
				bllin=c(bllin,rep(0,nrow(conp)))
				bulin=c(bllin,rep(0,nrow(conp)))
			}
		}
		for(i in 1:nparstotal) {if(fixed[i]==0) A[,i]=0}
		# check for zero rows in A which may result from fixing parameters, and remove them
		idx=numeric(0) 
		if(nrow(A)>0) {
			for(i in 1:nrow(A)) {if(sum(as.logical(A[i,]))!=0) idx = c(idx,i) }
			A=matrix(A[idx,],ncol=nparstotal)
		}
		#  ... also remove corresponding entries in bllin and bulin
		bllin=bllin[idx]
		bulin=bulin[idx]
		freepars=sum(as.logical(fixed))-qr(A[which(bllin==bulin),])$rank
		if(is.null(modname)) modname=paste(ng,"group model")
		mgd=list(ng=ng,modname=modname,mixmod=mixmod,nrcomp=nrcomp,nstates=nstates,itemtypes=itemtypes,
			npars=npars,nparstotal=nparstotal,freepars=freepars,
			freeparsnotd=freeparsnotd,pars=pars,fixed=fixed,
			A=A,bl=bl,bu=bu,bllin=bllin,bulin=bulin,td=td,tdfit=1,st=mixmod[[1]]$st)
	} #end ng>1
	if(nrcomp>1) class(mgd)=c("mgd","mixdmm")
	else class(mgd)=c("mgd","dmm")
	mgd
}

summary.mgd <- function(object, specs=FALSE, precision=3, se=NULL, ...) {
	if(object$ng==1) {
		if(class(object)[2]=="mixdmm") summary.mixdmm(object,specs,precision,se,...)
		else summary.dmm(object,specs,precision,se,...)
	} else {
		cat(" Model: ", object$modname, "\n")
		cat(" Nr of groups:     ", object$ng, "\n")
		cat(" Nr of parameters: ", object$nparstotal, "\n")
		cat(" Free parameters:  ", ifelse(object$tdfit,object$freepars,object$freeparsnotd), "\n")
		for(i in 1:object$ng) {
			cat(" Model for group: ", i, "\n")
			if(!is.null(se)) summary(object$mixmod[[i]],specs, precision, se=se[((i-1)*object$mixmod[[i]]$npars+1):((i)*object$mixmod[[i]]$npars)], ...)
			else summary(object$mixmod[[i]],specs, precision, ...)
		}
	}
}

