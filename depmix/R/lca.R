
# 
# LATENT CLASS MODEL (DMM WITHOUT TRANSITIONS)
# 

lca <- function(nclasses, itemtypes, modname=NULL, fixed=NULL, stval=NULL, conrows=NULL, conpat=NULL, linmat=NULL, snames=NULL, inames=NULL) {
# recode itemtypes/assign itemnames
	TYPES=c("gaussian","normal")
	itemtnames = sapply(itemtypes,FUN=function(x){pmatch(tolower(as.character(x)),TYPES)})
	itemtypes[which(!is.na(itemtnames))]=TYPES[itemtnames[which(!is.na(itemtnames))]]
	itemtnames=itemtypes
	itemtypes=recitt(itemtypes) #recitt in depmix-internal
	if(any(is.na(itemtypes))) stop("Itemtypes incorrectly specified (e.g. ga is underdetermined).")
	nitems = length(itemtypes)
# parameter counting
	lobs=sapply(itemtypes,FUN=np)
	lcmnp=nclasses*sum(lobs)+nclasses+1
# add appropriate numbers of zeroes/ones to fixed/stval/conrows/conpat
	if(!is.null(fixed)) {
		if(!(is.null(conpat))) warning("fixed will be overridden by conpat.")
		if(!(length(fixed)==lcmnp)) stop("fixed has incorrect length.")
		else fixed=c(0,rep(0,nclasses*nclasses),fixed[2:lcmnp])
	} else {
		fixed=c(0,rep(0,nclasses*nclasses),rep(1,(lcmnp-1)))
	}
	if(!is.null(stval)) {
		if(!(length(stval)==lcmnp)) stop("stval has incorrect length.")
		else stval=c(1,diag(nclasses),stval[2:lcmnp])
		st=TRUE
	} else {
		st=FALSE
		stval=c(1,diag(nclasses),runif((lcmnp-1)))
	}
	if(!is.null(conrows)) {
		conrows=matrix(conrows,ncol=lcmnp,byrow=TRUE)
		trzeroes=matrix(0,nrow=nrow(conrows),ncol=nclasses*nclasses)
		mpzeroes=matrix(0,nrow(conrows),1)
		conrows=t(cbind(mpzeroes,trzeroes,conrows[,2:(ncol(conrows))]))
	}
	if(!is.null(linmat)) {
		trzeroes=matrix(0,nrow=nrow(linmat),ncol=nclasses*nclasses)
		linmat=cbind(trzeroes,linmat)
	}
	if(!is.null(conpat)) {
		conpat=c(rep(0,nclasses*nclasses),conpat)
		fixed=NULL
	}
	if(is.null(snames)||length(snames)!=nclasses) snames=paste("Class",1:nclasses,sep="")
	if(is.null(modname)) modname=paste(nclasses,"class model")
	model=dmm(nstates=nclasses,itemtypes=itemtypes,modname=modname,fixed=fixed,stval=stval,conrows=conrows,conpat=conpat,linmat=linmat,snames=snames,inames=inames)
	class(model) = c("dmm","lcm")
	model$st=st
	model
}

# note that this is just a wrapper for dmm and it returns a full model with
# the transition parameters in there, however, they are simply discarded in
# estimation of the model, but they need to be considered when specifying 
# constraints