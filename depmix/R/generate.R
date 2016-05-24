
#####################################################
#                                             		#
#   DATA GENERATION ACCORDING TO (MIX)DMM OBJECT   	#
#                                             		#
#####################################################

generate <- function(ntimes,dmm,nreal=1) {
	if(class(dmm)[1]=="mgd") stop("run group models separately")
	ntimes = rep(ntimes,nreal)
	trans=list(); obser=list(); init=list()
	if(class(dmm)[1]=="dmm") dmm=mixdmm(dmm=list(dmm))
	for(i in 1:dmm$nrcomp) {
		idx=paridx(dmm$nstates,dmm$itemtypes,mat="tr",comp=i)
		trans[[i]]=matrix(dmm$pars[idx],dmm$nstates[i],byrow=TRUE)
		idx=paridx(dmm$nstates,dmm$itemtypes,mat="ob",comp=i)
		obser[[i]]=matrix(dmm$pars[idx],dmm$nstates[i],byrow=TRUE)
		idx=paridx(dmm$nstates,dmm$itemtypes,mat="in",comp=i)
		init[[i]]=matrix(dmm$pars[idx],dmm$nstates[i])
	}
	ntcount=0
	# create return value vectors
	nitems=length(dmm$itemtypes)
	obs = matrix(0,ncol=nitems,nrow=sum(ntimes))
	instates=integer(length(ntimes))
	incomp=integer(length(ntimes))
	for(j in 1:length(ntimes)) {
		# choose component
		if(dmm$nrcomp==1) ic=1
		else ic = which(rmultinom(1,size=1,prob=dmm$pars[1:dmm$nrcomp])==1)
		incomp[j]=ic
		# choose initial state
		ist = which(rmultinom(1,size=1,prob=init[[ic]])==1)
		st = ist
		instates[j]=ist
		for (t in 1:ntimes[j]) {
			obcount=0
			# choose observation values
			for(i in 1:nitems) {
				pars = obser[[ic]][st,(obcount+1):(obcount+np(dmm$itemtypes[i]))]
				obs[ntcount+t,i]=fresp(dmm$itemtypes[i],pars)
				obcount = obcount + np(dmm$itemtypes[i])
			}
			# choose transition
			nst = which(rmultinom(1,size=1,prob=trans[[ic]][st,])==1)
			st = nst
		}
		ntcount = ntcount + ntimes[j]
	}
	# repeat ntimes
	gen <- markovdata(dat=obs,itemtypes=dmm$itemtnames,ntimes=ntimes)
	attr(gen,"instates")=instates
	if(dmm$nrcomp>1) attr(gen,"incomp")=incomp
	gen
}
