# 
# DEPMIX INTERNAL FUNCTIONS, NOT TO BE CALLED BY USER
# 

# 
# AUXILIARIES FOR FITDMM
# 

checkSetRecode <- function(dat,dmm,tdcov=0,printlevel=1) {
	
## create the model to be fitted, changed into class mgd
	if(class(dmm)[1]=="dmm") {
 		xgmod <- mixdmm(dmm=list(dmm))	
		xgmod <- mgdmm(dmm=xgmod)
	}
	
	if(class(dmm)[1]=="mixdmm") xgmod <- mgdmm(dmm=dmm)
	if(class(dmm)[1]=="mgd") xgmod <- dmm
	

## some correctness checks, far from complete
	if(xgmod$ng==1) {
		if(!(class(dat)=="md")) stop("Data should have class md")
		dat <- list(dat)
	}
	
	if(length(dat)!=xgmod$ng) stop("Number of data sets does not match nr of groups.")
	if(!(length(xgmod$nstates)==xgmod$nrcomp)) stop("Length of nstates and nrcomp do not match.")
	if(!(length(xgmod$itemtypes)==(length(itemtypes(dat[[1]]))-ncov(dat[[1]])))) stop("Model incompatible with data, different numbers of items.")
	
	ncov=ncov(dat[[1]])
	if(tdcov==1 && ncov==0) stop("There are no covariates in the data so these cannot be used.")
	if(tdcov) if(xgmod$td==FALSE) stop("Model does not specify covariates.")
	
	if(printlevel>19) print("All checks okay.")
	
	gc()
	
## split data and covariates, data may have covariates, even when they're not fitted,
## so the data allways has to be split
	if(ncov>0) {
		nrit=dim(dat[[1]])[2]-ncov
		itt=rep(1,ncov)
		dcov=list()
		datsplit=list()
		for(i in 1:xgmod$ng) {
			dcov[[i]]=markovdata(dat=dat[[i]][,(nrit+1):(dim(dat[[i]])[2]),drop=FALSE],itemtypes=itt,ntimes=ntimes(dat[[i]]))
			datsplit[[i]]=markovdata(dat=dat[[i]][,1:nrit,drop=FALSE],itemtypes=itemtypes(dat[[i]])[1:nrit],ntimes=ntimes(dat[[i]]),replicates=replicates(dat[[i]]))
		}
		dat=datsplit
		if(printlevel>19) print("Data split okay.")
	}
	
## check missing data codes
	xm=-9999
	for(i in 1:xgmod$ng) {
		dat[[i]]=replace(dat[[i]],which(is.na(dat[[i]])),xm)
		if(ncov>0) datsplit[[i]]=replace(datsplit[[i]],which(is.na(datsplit[[i]])),xm)
	}
	
	if(printlevel>19) print("Missing data codes set.")
	
## send td model to C routines	
	if(tdcov) {
		nitems=length(xgmod$itemtypes)
 		z <- .C("covSetUp",
 			as.integer(xgmod$ng),
 			as.integer(xgmod$nrcomp),
 			as.integer(xgmod$nstates),
 			as.integer(nitems),
 			as.integer(xgmod$itemtypes),
 			as.double(xgmod$pars[(xgmod$npars+1):(2*xgmod$npars)]),
 			as.integer(xm),
 			as.integer(printlevel),
 			PACKAGE="depmix")

		## ... and set covariates in C
		mcov <- .C("multiCovSetUp",
			as.integer(length(dcov)),
			as.integer(printlevel),
			PACKAGE="depmix")
		datset <- list()
		for(i in 1:xgmod$ng) {
# 			catit=0
# 			dc=numeric(0)
			datset[[i]] <- .C("ngCovSetUp",
				as.integer(i),
				as.double(t(dcov[[i]])),
				as.integer(dim(dcov[[i]])[2]),
				as.integer(itemtypes(dcov[[i]])),
				as.integer(length(ntimes(dcov[[i]]))),
				as.integer(ntimes(dcov[[i]])),
				as.double(replicates(dcov[[i]])),
				as.integer(xm),
				as.integer(printlevel),
				PACKAGE="depmix")
		}
		if(printlevel>19) print("Td model set in C")
	}
	
 	nitems=length(xgmod$itemtypes)
 	pars=xgmod$pars[1:xgmod$npars] # only the normal model pars, tdpars are already set
	
	if(printlevel>19) print("Preparing to send model to C routines.")
	
## send model to C routines
	z <- .C("mixModelSetUp",
		as.integer(xgmod$ng),
		as.integer(xgmod$nrcomp),
		as.integer(xgmod$nstates),
		as.integer(nitems),
		as.integer(xgmod$itemtypes),
		as.double(pars),
		as.integer(xm),
		as.integer(printlevel),
		PACKAGE="depmix")
	
	if(printlevel>19) print("Model set up finished in C.")
		
## send data to C routines
	mdat <- .C("multiDataSetUp",
		as.integer(length(dat)),
		as.integer(printlevel),
		PACKAGE="depmix")
	datset <- list()
	
	for(i in 1:xgmod$ng) {
  		attr(dat[[i]],"itemtypes")=replace(itemtypes(dat[[i]]),which(itemtypes(dat[[i]])!="categorical"),1)
		itt=rep(0,dim(dat[[i]])[2])
		## recode itemtypes 
		for(it in 1:dim(dat[[i]])[2]) {
			if(itemtypes(dat[[i]])[it]=="categorical") itt[it]=length(setdiff(unique(dat[[i]][,it]),xm))	
			else itt[it]=1
		}
 		attr(dat[[i]],"itemtypes")=itt
		# recode the data if necessary
		catit=0
		dc=numeric(0)
		for(it in 1:dim(dat[[i]])[2]) {
			if(itemtypes(dat[[i]])[it]>1) {
 				catit=catit+1
 				datacats=sort(setdiff(unique(dat[[i]][,it]),xm))
 				if(all.equal(datacats,(1:itt[it]))==TRUE) dc[catit]=1
 				else dc[catit]=0
			}
		}
		if(all(as.logical(dc))) recdat=dat[[i]]
		else recdat=recode(dat[[i]],xm)
		
		datset[[i]] <- .C("ngDataSetUp",
			as.integer(i),
			as.double(t(recdat)),
			as.integer(dim(dat[[i]])[2]),
			as.integer(itemtypes(dat[[i]])),
			as.integer(length(ntimes(dat[[i]]))),
			as.integer(ntimes(dat[[i]])),
 			as.double(replicates(dat[[i]])),
			as.integer(xm),
			as.integer(printlevel),
			PACKAGE="depmix")
	}
	
	if(printlevel>19) print("Data set up finished in C\n")
	
	return(xgmod)
}

recode <- function(dat,xm) {
	recdat=dat
	nit=length(itemtypes(dat))
	for(it in 1:nit) {
		if(itemtypes(dat)[it]>1) {
			original=sort(setdiff(unique(dat[,it]),xm))
			normvals=1:itemtypes(dat)[it]
			for(i in 1:itemtypes(dat)[it]) {
				recdat[which(dat[,it]==original[i]),it]=normvals[i]
			}
		}
	}
	recdat
}

# 
# AUXILIARIES FOR DMM
# 

# function fblo gives lower bounds of parameters for each density
fblo<-function(x,i,bigB) {
	if(x > 1)  return(rep(0,x))				# multinomial
	if(x== 1)  return(c(-bigB,-bigB))		# normal (mu - mean, sigma - stddev)
}

# function fbuo gives upper bounds of parameters for each density
fbuo<-function(x,i,bigB) {
	ubshift=+bigB	
	if(x > 1)  return(rep(1,x)) 				# multinomial
	if(x==1)   return(c(+bigB,+bigB)) 			# normal (mu - mean, sigma - stddev)
}

# function ppar generates random start values of parameters for each density and pastes stvals together
ppar<-function(x,z) {
	if(x > 1)  {y<-runif(x);y<-y/sum(y)} 							# multinomial
	if(x==1)    y<-c(rnorm(1), rnorm(1,1,0.1)) 						# normal (mu - mean, sigma - stddev)
	return(c(z,y)) 
}

#recode itemtypes from characters to magic numbers
recitt <- function(itemtypes) {
	itemtypes[itemtypes=="gaussian"|itemtypes=="normal"|itemtypes=="norm"]=1
	itemtypes=as.numeric(itemtypes)
	itemtypes
}

# function pp generates vector with parameter names for each itemtype
pp<-function(x){
	if(x > 1)  return(paste("p",1:x)) # multinomial
	if(x== 1)  return(c("mean", "stddev"))  # normal
}

# returns the number of parameters for each density
np <- function(x) {
	if(x > 1)  {
		return(x) 
	} else {
		if(x==1) {
			return(2)
		} else {
			stop("Unknown itemtype code in function np.")
		}
	}
}

# convert conpat vector to rows of constraint matrix
pa2conr <- function(x) {
	fix=as.logical(x)
	x=replace(x,list=which(x==1),0)
	un=setdiff(unique(x),0)
	y=matrix(0,0,length(x))
	for(i in un) {
		z=which(x==i)
		for(j in 2:length(z)) {
			k=rep(0,length(x))
			k[z[1]]=1
			k[z[j]]=-1
			y=rbind(y,k)
		}
	}
	pa = list(fix=fix,conr=y)
	return(pa)
}

##################################
#                                #
#   PARAMETER INDEX GENERATION   #
#                                #
##################################

paridx <- function(nstates,itemtypes,mat,idx1=0,idx2=0,it=0,comp=1,group=1) {	
#   all denotes all parameters of a component model, ie not for a group, and hence
# 	it does not include the mixing proportions
	MATS=c("mix","trans","obser","init","all") 
	mat = pmatch(mat,MATS)
	nitems=length(itemtypes)
	nrcomp=length(nstates)
	nstates=as.vector(nstates)
	## nr of pars per item, the sum of this is the nr of pars per state 
	# for the response parameters
	lobs=sapply(itemtypes,FUN=np)
	## nr of pars per component
	nrcp = rep(0,nrcomp)
	for(i in 1:nrcomp) nrcp[i]=nstates[i]*(nstates[i]+1)+sum(lobs)*nstates[i]
	ngpars=nrcomp+sum(nrcp)
	mix=1:nrcomp
	# mix prop pars
	if(mat==1) if(idx1==0) idx=mix else idx=mix[idx1]
	# trans pars
	if(mat==2) {
		tridx <- matrix(1:(nstates[comp]*nstates[comp]),nstates[comp],nstates[comp],byrow=TRUE)
		if(idx1==0) idx=c(t(tridx))
		else { if(idx2==0) { 
				idx=tridx[idx1,]
			}
			else idx=tridx[idx1,idx2]
		}
	}
	# obser pars
	if(mat==3) {
		en=nstates[comp]*nstates[comp]
		obidx=matrix((en+1):(nstates[comp]*sum(lobs)+en),nstates[comp],sum(lobs),byrow=TRUE)
		if(it==0) {
			if(idx1==0) { idx=c(t(obidx)) }
			else {
				if(idx2==0) { 
					idx=obidx[idx1,] 
				}
				else idx=obidx[idx1,idx2]
			}
		} else { # it>0
			if(nitems>1 && it>1) npprev=sum(sapply(itemtypes[1:(it-1)],FUN=np))
			else npprev=0
			if(itemtypes[it]==1) np=2 
			else np=itemtypes[it]
			if(idx1==0) { #return indices of item parameters for all states
				idx=c(t(obidx[,(npprev+1):(npprev+np)]))
			} else {
				if(idx2==0) idx=obidx[idx1,(npprev+1):(npprev+np)]
				else idx=obidx[idx1,npprev+idx2]
			}
		}
	}
	# init pars
	if(mat==4) {
		inidx=(nrcp[comp]-nstates[comp]+1):(nrcp[comp])
		if(idx1==0) { idx=inidx }
		else idx=inidx[idx1]
	}
	if(mat==5) idx=1:nrcp[comp]
 	if(comp>1) idx=idx+sum(nrcp[1:(comp-1)])
	if(mat!=1) idx=idx+nrcomp
	if(group>1) idx = idx + (group-1)*ngpars
	return(idx)
}

# function fresp generates random deviates for itemtype x, given 'pars', ie the item distribution parameters
fresp <- function(x,pars) {
	if(x > 1)  resp = which(rmultinom(1,size=1,prob=pars)==1)	# multinomial
	if(x==1)   resp = rnorm(1,pars[1],pars[2]) 					# normal (mu - mean, sigma - stddev)
	return(resp) 
}

# convert a list of matrices x to a block diagonal matrix 
# shamelessly stolen from package assist
bdiag <- function (x) {
    if (is.matrix(x)) 
        return(x)
    if (!is.list(x)) 
        stop("dismatched input")
    n <- length(x)
    x <- lapply(x, function(y) if (length(y)) 
        as.matrix(y)
    else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1, ]
    cc <- d[2, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1, -1] <- rcum[-n]
    ind[2, ] <- rcum
    ind[3, -1] <- ccum[-n]
    ind[4, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) {
        imat[(y[1] + 1):y[2], (y[3] + 1):y[4]]
    }, imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    out
}

# 
# STARTING VALUE GENERATION (AUXILIARY) FUNCTIONS
# 

kmstart <- function(dat,dmm) {
	nst=dmm$nstates
	kmc=try(kmeans(na.omit(dat[,1:(nitems(dat)-ncov(dat))]),nst)$cluster,silent=TRUE)
	if(class(kmc)=="try-error") {
		while(class(kmc)=="try-error") {
			kmc=try(kmeans(dat[,1:(nitems(dat)-ncov(dat))],nst)$cluster,silent=TRUE)
		}	
	}
	st=cl2st(kmc,dat,dmm)
	st
}

poststart <- function(dat,dmm) {
	post=posterior(dat,dmm)
	st=cl2st(post$states[[1]][,2],dat,dmm)
	st
}

cl2st <- function(cluster,dat,dmm) {
	nst=dmm$nstates
	stin=numeric(nst)
	if("lcm"%in%class(dmm) || dmm$nst==1) {
		sttr=diag(nst)
		for(i in 1:nst) {
			stin[i]=length(cluster[which(cluster==i)]) 
		}
	} else {
		sttr=matrix(0,nst,nst)
		for(i in 1:ind(dat)) {
			bg=ifelse(i>1,sum(ntimes(dat)[1:(i-1)])+1,1)
			en=bg+ntimes(dat)[i]-1
			# add diagonal elements to make sure that the transition table is complete, which
			# is not neccessarily so for each case, ie not all transitions necc occur
			cl1=c(as.numeric(gl(nst,2)),cluster[bg:en])
			cl2=cl1[-1]
			cl1=cl1[-length(cl1)]
   			sttr=sttr+matrix(table(cl1,cl2),nst,nst)
			stin[cluster[bg]]=stin[cluster[bg]]+1
			# correct for the added transition from nst to the start state
   			sttr[nst,cluster[bg]]=sttr[nst,cluster[bg]]-1
		}
		# correct for the added transitions between i and i+1
   		for(i in 1:nst) {
			sttr[i,i]=sttr[i,i]-ind(dat)
			if(i<nst) sttr[i,i+1]=sttr[i,i+1]-ind(dat)
		}
		sttr=sttr/rowSums(sttr)
	}
	if(any(stin==0)) stin=stin+0.05
	stin=stin/sum(stin)
	stob=cl2stob(cluster,dat,dmm)
	st=c(1,t(sttr),stob,stin)
	st
}

tr2stin <- function(sttr) {
	etr=eigen(t(sttr))
	stin=etr$vectors[,which(etr$values==1)]
	stin=stin/sum(stin)
	stin
}
	
cl2stob <- function(cluster,dat,dmm) {
	stob=numeric(0)
	nst=dmm$nstates
	for(i in 1:nst) {
		for(j in 1:length(dmm$itemtypes)) {
			if(dmm$itemtypes[j]==1) {
				meanob=mean(dat[which(cluster==i),j])
				sdob=sd(dat[which(cluster==i),j])
				stob=c(stob,meanob,sdob)
			}
			if(dmm$itemtypes[j]>1) {
				y=sort(unique(na.omit(dat[,j])))
				catprobs=numeric(dmm$itemtypes[j])
				for(k in 1:length(y)) {
					catprobs[k]=length(dat[which(cluster==i & dat[,j]==y[k]),j])
				}
				catprobs=catprobs/sum(catprobs)
				stob=c(stob,catprobs)
			}
		}
	}
	stob
}