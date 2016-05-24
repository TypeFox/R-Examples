#################################################
#                                               #
#   DEPENDENT MIXTURE MODEL OBJECT DEFINITION   #
#                                               #
#################################################

dmm <- function(nstates, itemtypes, modname=NULL, fixed=NULL, stval=NULL, conrows=NULL, conpat=NULL, tdfix=NULL, tdst=NULL, linmat=NULL, snames=NULL, inames=NULL) { 
# check on nstates
	if(nstates<1) stop("nstates must be positive.")
# recode itemtypes/assign itemnames
	TYPES=c("gaussian","normal")
	itemtnames = sapply(itemtypes,FUN=function(x){pmatch(tolower(as.character(x)),TYPES)})
	itemtypes[which(!is.na(itemtnames))]=TYPES[itemtnames[which(!is.na(itemtnames))]]
	itemtnames=itemtypes
	itemtypes=recitt(itemtypes) #recitt in depmix-internal
	if(any(is.na(itemtypes))) stop("Itemtypes incorrectly specified (e.g. ga is underdetermined).")
	nitems = length(itemtypes)
# parameter counting
	# nr of pars per item, the sum of this is the nr of pars per state for the response parameters
	lobs=sapply(itemtypes,FUN=np) # np in depmix-internal
	npars=nstates*nstates+nstates*sum(lobs)+nstates+1
# checks on fixed values, start values and conpat
	fv=FALSE
	if(!is.null(fixed)) fv=TRUE
	if(!is.null(conpat)) {
		if(fv) warning("Fixed will be overridden by conpat, see ?dmm.")
		if(!(length(conpat)==npars)) stop("Conpat has incorrect length, see ?dmm.")
		fixed=pa2conr(conpat)$fix
		conp=pa2conr(conpat)$conr #pa2conr in depmix-internal
		fv=TRUE
		if(nrow(conp)>0) pa=TRUE
		else pa=FALSE
	}
	else pa=FALSE
	if(is.null(stval)) {
		if(fv) stop("When there are fixed parameters, start values should be provided in stval.")
		st=FALSE
	} else {
		if(length(stval)!=npars) stop("Stval has incorrect length (it should include the mixing proportion=1.0), see ?dmm.")
		else st=TRUE
	}
	cr=FALSE
	if(!is.null(conrows)) {
		if((length(conrows)%%npars)!=0) stop("Length of conrows is not a multiple of the number of parameters.") 
		else cr=TRUE
	}
	if(!fv) fixed=c(0,rep(1,npars-1))
 	if(fixed[1]!=0) {fixed[1]=0; warning("The mixing proportion of a single chain is fixed to 1.")}
	if(nstates==1) {fixed[2]=0; fixed[npars]=0}
	bigB <- 10^21 #upper and lower limit infinite, npsol magical number
# provide parameter start values if none are given, scale start values if provided
	mixprop=1.0 # mixing proportion (added for consistency)
	if(st) {
		trans <- matrix(stval[2:(nstates*nstates+1)],nrow=nstates,byrow=TRUE)
		obser <- matrix(stval[(nstates*nstates+2):(npars-nstates)],nrow=nstates,byrow=TRUE)
		init=stval[(npars-nstates+1):npars]
	} else {
		trans <- matrix(runif(nstates*nstates,0,1),nrow=nstates)
		obser <- matrix(0,nrow=nstates,ncol=sum(lobs))
		for(j in 1:nstates) {
			z=numeric(0)
			for (i in 1:nitems) z=ppar(itemtypes[i],z)
			obser[j,]=z
		} #ppar in depmix-internal
		init=runif(nstates,0,1)
	}
# scale init & trans parameters parameters
	trans <- trans/rowSums(trans) 
	init=init/sum(init)
# join parameters together
	pars=c(mixprop,t(trans),t(obser),init)
# scale observation parameters for itemtypes>1
	for(j in 1:nstates) {
		for(i in 1:nitems) {
			if(itemtypes[i]>1) { 
 				itpars=pars[paridx(nstates,itemtypes,mat="ob",idx1=j,it=i)]
				pars[paridx(nstates,itemtypes,mat="ob",idx1=j,it=i)] = itpars/sum(itpars)
			}
		}
	}
# bounds for transition pars
	blt <- rep(0,nstates*nstates)
	but <- rep(1,nstates*nstates)
# bounds for obser pars
	blo <- numeric(0); buo <- numeric(0)
	for(j in 1:nstates) { 
		for(i in 1:nitems) {
			blo=c(blo,fblo(itemtypes[i],i,bigB))
			buo=c(buo,fbuo(itemtypes[i],i,bigB))
		}
	}
# bounds for init parameters
	bli=rep(0,nstates)
	bui=rep(1,nstates)
# make bl and bu (first 1's for the mixing proportion)
	bl <- c(1.0,blt,blo,bli) 
	bu <- c(1.0,but,buo,bui)
# set bounds for fixed parameters
	bl[which(fixed==0)]=pars[which(fixed==0)]
	bu[which(fixed==0)]=pars[which(fixed==0)]
# trans matrix sum constraints
	blst=numeric(0); bust=numeric(0)
	if(nstates>1) {
		sctr=matrix(0,nstates,npars)
		for(i in 1:nstates) sctr[i,paridx(nstates,itemtypes,mat="tr",i)]=1
		blst <- rep(1,nstates)
		bust <- rep(1,nstates)
	}
# obser matrix sum constraints
	scob=matrix(0,0,npars); blso=numeric(0); buso=numeric(0)
	for(i in 1: nstates) {
		for(j in 1:nitems) {
			if(itemtypes[j]>1) {
				sc=rep(0,npars)
				sc[paridx(nstates,itemtypes,mat="ob",it=j,idx1=i)]=1
				scob=rbind(scob,sc)
				blso=c(blso,1); buso=c(buso,1)
			}
		}
	}
# init vector sum constraint
	blsi=numeric(0); busi=numeric(0)
	if(nstates>1) {
		scin=rep(0,npars)
		scin[(npars-nstates+1):(npars)]=1
		blsi=c(1)
		busi=c(1)
	}
# join linear constraint bounds and make linear constraint matrix A
	bllin <- c(blst,blso,blsi)
	bulin <- c(bust,buso,busi)
	if(nstates==1) A=scob
	else A=rbind(sctr,scob,scin)
# add user supplied linear constraints to A and set their bounds
	if(cr) {
		conrows <- matrix(conrows,ncol=npars,byrow=TRUE)
		A=rbind(A,conrows)
		bllin <- c(bllin,rep(0,nrow(conrows)))
		bulin <- c(bulin,rep(0,nrow(conrows)))
	}
	if(pa) {
		A=rbind(A,conp)
		bllin <- c(bllin,rep(0,nrow(conp)))
		bulin <- c(bulin,rep(0,nrow(conp)))
	}
 	freeparsnotd=sum(as.logical(fixed))
 	Ared=A[,which(fixed==1),drop=FALSE]
 	freeparsnotd=freeparsnotd-qr(Ared[which(bllin==bulin),,drop=FALSE])$rank
# time-dependent covariates
	td=0; tdtr=FALSE; tdob=FALSE; tdin=FALSE
	if(!is.null(tdfix)) {
		td=1 # only one covariate supported at the moment 
		bltd=rep(0,npars)
		butd=rep(0,npars)
		if(any(tdfix[paridx(nstates,itemtypes,mat="tr")]==1)) {
			tdtr=TRUE
			if(nstates==1) stop("Covariates on the transition matrix do not make sense with a single state model.\n")
			butd[paridx(nstates,itemtypes,mat="tr")]=1
			bltd[paridx(nstates,itemtypes,mat="tr")]=-1
		}
		# etc for other betas
		if(any(tdfix[paridx(nstates,itemtypes,mat="ob")]==1)) {
			tdob=TRUE
			butd[paridx(nstates,itemtypes,mat="ob")]=bigB
			bltd[paridx(nstates,itemtypes,mat="ob")]=-bigB
		}
		if(any(tdfix[paridx(nstates,itemtypes,mat="in")]==1)) {
			tdin=TRUE
			if(nstates==1) stop("Covariates on the initial probs do not make sense with a single state model.\n")
			butd[paridx(nstates,itemtypes,mat="in")]=1
			bltd[paridx(nstates,itemtypes,mat="in")]=-1
		}
		tdpars=rep(0,npars)
		if(is.null(tdst)) tdpars=replace(tdpars,tdfix==1,rnorm(sum(tdfix),0,0.1))
		else tdpars=tdst
	}
# covariates for tr parameters
	if(tdtr) {
		# values: make values conformable with constraints
		tdtrans=matrix(tdpars[paridx(nstates,itemtypes,mat="tr")],nrow=nstates,byrow=TRUE)
		if(nstates==2) tdtrans[,2]=-tdtrans[,1]
		else tdtrans[,nstates]=-apply(tdtrans[,(1:(nstates-1))],1,sum)
		tdpars[paridx(nstates,itemtypes,mat="tr")]=t(tdtrans)
		# make the constraints
		sctdtr=matrix(0,nstates+nstates*nstates,2*npars)
		blsctdtr=rep(0,nstates+nstates*nstates)
		busctdtr=c(rep(0,nstates),rep(1,nstates*nstates))
		# ... between the beta's (these sum to zero per state)
		for(i in 1:nstates) sctdtr[i,paridx(nstates,itemtypes,mat="tr",idx1=i)+npars]=1
		# ... between transpars and beta's  0<=a_ij+b_ij<=1, assuming the covariate is between 0 and 1
 		for(i in 1:nstates) {
 			for(j in 1:nstates) {
 				pidx=paridx(nstates,itemtypes,mat="tr",idx1=i,idx2=j)
				sctdtr[nstates+(i-1)*nstates+j,c(pidx,pidx+npars)]=c(1,1)
 			}
 		}
	}
# covariates for obser parameters
	if(tdob) {
		# this is the constraint matrix with enough rows for the zero sum constraints
		sctdob=matrix(0,nstates*(sum(as.logical(itemtypes>1))+sum(itemtypes[itemtypes>1])),2*npars)
		blsctdob=c(rep(0,nstates*(sum(as.logical(itemtypes>1)))),rep(0,(nstates*(sum(itemtypes[itemtypes>1])))))
		busctdob=c(rep(0,nstates*(sum(as.logical(itemtypes>1)))),rep(1,(nstates*(sum(itemtypes[itemtypes>1])))))
		# values are already given above, but they still have to be made consistent with the constraints
		kk=0
		for(i in 1:nstates) { 
			for(j in 1:nitems) {
				if(itemtypes[j]>1) {
					kk=kk+1
					# these beta's have to sum to one
					pp=tdpars[paridx(nstates,itemtypes,mat="ob",idx1=i,it=j)]
					pp[itemtypes[j]]=-sum(pp[1:(itemtypes[j]-1)])
					tdpars[paridx(nstates,itemtypes,mat="ob",idx1=i,it=j)]=pp
					sctdob[kk,paridx(nstates,itemtypes,mat="ob",idx1=i,it=j)+npars]=1
				}
			}
		}
		# constraints
		# additionally obser_ij + beta_ij have to sum to one together for each state i and category j
		for(i in 1:nstates) { 
			for(j in 1:nitems) {
				if(itemtypes[j]>1) {
					for(k in 1:itemtypes[j]) {
						kk=kk+1
						pidx=paridx(nstates,itemtypes,mat="ob",idx1=i,it=j,idx2=k)
						sctdob[kk,c(pidx,pidx+npars)]=c(1,1)
					}
				}
			}
		}
	}
# covariates for init paramaters
	if(tdin) {
		tdinit=tdpars[(npars-nstates+1):npars]
		# values have to sum to zero
		if(nstates==2) tdinit[2]=-tdinit[1]
		else tdinit[nstates]=-sum(tdtrans[1:(nstates-1)])
		# make constraint matrix rows
		sctdin=matrix(0,1+nstates,2*npars)
		blsctdin=rep(0,1+nstates)
		busctdin=c(0,rep(1,nstates))
		# beta's sum to zero
		sctdin[,paridx(nstates,itemtypes,mat="in")+npars]=1
		# beta + initpar should be between zero and one
		for(i in 1:nstates) {
			pidx=paridx(nstates,itemtypes,mat="in",idx1=i)
			sctdin[1+i,c(pidx,pidx+npars)]=c(1,1)
		}
	}
# add contraints for the covariate parameters
	nparstotal=npars
	if(td) {
		nparstotal=2*npars
		pars=c(pars,tdpars)
		fixed=c(fixed,tdfix)
		bu=c(bu,butd)
		bl=c(bl,bltd)
		betaA=matrix(0,nrow(A),npars)
		A=cbind(A,betaA)
		if(tdtr) {
			bllin=c(bllin,blsctdtr)
			bulin=c(bulin,busctdtr)
			A=rbind(A,sctdtr)
		}
		if(tdob) {
			bllin=c(bllin,blsctdob)
			bulin=c(bulin,busctdob)
			A=rbind(A,sctdob)
		}
		if(tdin) {
			bllin=c(bllin,blsctdin)
			bulin=c(bulin,busctdin)
			A=rbind(A,sctdin)
		}
	}
	# set bounds for fixed pars
	if(nparstotal>npars) {
		bl[which(fixed==0)]=pars[which(fixed==0)]
		bu[which(fixed==0)]=pars[which(fixed==0)]
	}
	# check for zero rows in A which may result from fixing parameters at zero, and remove them
	# set zeroes in A when pars are fixed, and change bllin and bulin accordingly
	if(nrow(A)>0) {
		Bcomplete=matrix(as.logical(A),nrow=nrow(A))
		xcomplete=apply(Bcomplete,1,sum)
		A[,which(fixed==0)]=0
		B=matrix(as.logical(A),nrow=nrow(A))
		x=apply(B,1,sum)
		for(i in 1:nrow(A)) {
			if(x[i]!=xcomplete[i]) {
				for(j in 1:nparstotal) {
					if(fixed[j]==0 && Bcomplete[i,j]!=0 && pars[j]!=0) {
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
	if(!is.null(linmat)) A <- linmat
	# count free parameters
 	freepars=sum(as.logical(fixed))
	Ared=A[,which(fixed==1),drop=FALSE]
	freepars=freepars-qr(Ared[which(bllin==bulin),,drop=FALSE])$rank
# names
	if(is.null(modname)) modname=paste(nstates,"-state model")
	if(is.null(snames)||length(snames)!=nstates) snames=paste("State",1:nstates,sep="")
	if(is.null(inames)||length(inames)!=nitems) inames=rep(paste("Item",1:nitems,sep=""),lobs)  
# model (return value)
	model=list(modname=modname,nstates=nstates,snames=snames,nitems=nitems,
	itemtypes=itemtypes,itemtnames=itemtnames,inames=inames,
	npars=npars,nparstotal=nparstotal,freepars=freepars,
	freeparsnotd=freeparsnotd,pars=pars,fixed=fixed,
	A=A,bl=bl,bu=bu,bllin=bllin,bulin=bulin,
	td=td,tdtr=tdtr,tdob=tdob,tdin=tdin,tdfit=1,st=st) 
	class(model) = "dmm"
	model
}

summary.dmm <- function(object, specs=FALSE, precision=3, se=NULL, ...) {
	cat(" Model: ", object$modname, " \n")
	cat(" Number of parameters: ", object$nparstotal, " \n")
	cat(" Free parameters:      ", ifelse(object$tdfit,object$freepars,object$freeparsnotd), " \n")
	cat(" Number of states:     ", object$nstates,"\n")
	cat(" Number of items:      ", object$nitems,"\n")
	cat(" Item types:           ", object$itemtnames,"\n\n")
	ses=ifelse(is.null(se),0,1)
	if(!"lcm" %in% class(object)) {
		cat(" Parameter values, transition matrix \n\n")
		trsx=matrix(object$pars[paridx(object$nstates,object$itemtypes,mat="tr")],object$nstates,byrow=TRUE)
		rownames(trsx)=object$snames
		if(object$tdtr && object$tdfit==1) {
            betr=matrix(object$pars[paridx(object$nstates,object$itemtypes,mat="tr")+object$npars],object$nstates,byrow=TRUE)
            trsx=matrix(rbind(c(trsx),c(betr)),object$nstates*2,object$nstates)
            be=rep(x="be", object$nstates)
            rownames(trsx)=c(rbind(object$snames,be))
        }
		if(ses) {
			setr=matrix(se[paridx(object$nstates,object$itemtypes,mat="tr")],object$nstates,byrow=TRUE)
			tr=matrix(0,object$nstates*3,object$nstates)
			for(i in 1: object$nstates) {
				tr[3*(i-1)+1,]=trsx[i,]
				tr[3*(i-1)+2,]=setr[i,]
				for(j in 1: object$nstates) tr[3*(i-1)+3,j]=ifelse(tr[3*(i-1)+2,j]==0,NA,tr[3*(i-1)+1,j]/tr[3*(i-1)+2,j])
			}
			trsx=tr
			pn=rep(x="se", object$nstates)
			tn=rep(x="t", object$nstates)
			rownames(trsx)=c(rbind(object$snames,pn,tn))
		}
		colnames(trsx)=object$snames
		print(round(trsx,precision))
		cat("\n\n")
	}
	cat(" Parameter values, observation parameters \n\n")
	obsx=matrix(object$pars[paridx(object$nstates,object$itemtypes,mat="ob")],object$nstates,byrow=TRUE)
	rownames(obsx)=object$snames
	if(object$tdob && object$tdfit==1) {
		tdobpars=matrix(object$pars[paridx(object$nstates,object$itemtypes,mat="ob")+object$npars],object$nstates,byrow=TRUE)
		obsx=matrix(rbind(c(obsx),c(tdobpars)),object$nstates*2,byrow=FALSE)
		be=rep(x="be", object$nstates)
		rownames(obsx)=c(rbind(object$snames,be))
	}
	if(ses) {
		seob=matrix(se[paridx(object$nstates,object$itemtypes,mat="ob")],object$nstates,byrow=TRUE)
		ob=matrix(0,object$nstates*3,ncol=ncol(obsx))
		for(i in 1: object$nstates) {
			ob[3*(i-1)+1,]=obsx[i,]
			ob[3*(i-1)+2,]=seob[i,]
			for(j in 1: ncol(obsx)) ob[3*(i-1)+3,j]=ifelse(ob[3*(i-1)+2,j]==0,NA,ob[3*(i-1)+1,j]/ob[3*(i-1)+2,j])
		}
		obsx=ob
		pn=rep(x="se", object$nstates)
		tn=rep(x="t", object$nstates)
		rownames(obsx)=c(rbind(object$snames,pn,tn))
	}
	lobs=sapply(object$itemtypes,FUN=np)
	inames=object$inames
	pnames=unlist(sapply(object$itemtypes,FUN=pp))
	nam=paste(inames,pnames,sep=",")
	colnames(obsx)=nam
	print(round(obsx,precision))
	cat("\n\n")
	if(!"lcm" %in% class(object)) cat(" Parameter values, initial state probabilies \n\n")
	else cat(" Parameter values, unconditional (class) probabilities \n\n")
	if(ses) {
		inix=matrix(c(object$pars[paridx(object$nstates,object$itemtypes,mat="in")],rep(0,2*object$nstates)),3,object$nstates,byrow=TRUE)
		sein=se[paridx(object$nstates,object$itemtypes,mat="in")]
		for(i in 1: object$nstates) {
			inix[2,i]=sein[i]
			inix[3,i]=ifelse(sein[i]==0,NA,inix[1,i]/sein[i])
		}
		rownames(inix)=c("val","se","t")
	} else {
		inix=matrix(object$pars[paridx(object$nstates,object$itemtypes,mat="in")],1,object$nstates,byrow=TRUE)
		rownames(inix)="val"
 		if(object$tdin && object$tdfit==1) {
			inix=rbind(inix,object$pars[paridx(object$nstates,object$itemtypes,mat="in")+object$npars])
			rownames(inix)=c("val","be")
		}
	}
	colnames(inix)=object$snames
	print(round(inix,precision))
	cat("\n")
}
