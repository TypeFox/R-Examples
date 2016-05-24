"opt.par" <-
function(
	#function for optimizig partition in blockmodeling
	M,	#matrix
	clu,	#partition or a list of paritions for each mode
#	e1,	#weight of the error of binearized matrix
#	e2,	#weight of the error of valued conenctions
	#m=NULL,	#suficient value individual cells
#	s="default",	#suficient value for colum and row statistics
#	FUN,	#function to calculate row and colum statistics
#	blocks#=c("null","com","reg"),	#permissible block types and their ordering
#	BLOCKS=NULL,	#array of permissible block types and their ordering for all blocks
#	mindim = 2,	#minimal dimension for regulal, dominant and functional blocks
#	save.err.v=FALSE,	#save a vector of errors of all block types for all blocks
	approach,
	...,	#other arguments to called functions - to 'crit.fun'  or 'opt.par.tmp'
	maxiter=50,	#maximum number of iterations
	trace.iter=FALSE,	#save a result of each iteration or only the best (minimal error)
	switch.names=NULL,	#should partitions that only differ in group names be considert equal (is c(1,1,2)==c(2,2,1))
	save.initial.param=TRUE,	#should the initial parameters be saved
#	force.fun=NULL, #select the function used to evaluate partition
	skip.par=NULL,	#the partions that are not allowed or were already checked and should be skiped
	save.checked.par=!is.null(skip.par),	#should the checked partitions be saved
	merge.save.skip.par=all(!is.null(skip.par),save.checked.par), #should the checked partitions be merged with skiped ones
	check.skip="never"	#when should the check be preformed:
								# "all"  - before every call to 'crit.fun'
								# "iter" - at the end of eack iteratiton
								# "opt.par"  - before every call to 'opt.par', implemented in opt.these.par and opt.random.par
								# "never" - never
){
	if(is.list(clu)){
		k<-sapply(clu,function(x)length(unique(x)))
	} else k<-length(unique(clu))
	
	dots<-list(...)
	if(is.null(switch.names)){
		switch.names<-is.null(dots$BLOCKS)
	}
	
	
	nmode<-length(k)

	if(nmode==1) clu<-as.integer(factor(clu))
	if(nmode==2) clu<-lapply(clu,function(x)as.integer(factor(x)))
	if(nmode>2) {
		clu<-lapply(clu,function(x)as.integer(factor(x)))
		for(i in 2:length(clu)){
			clu[[i]]<-clu[[i]] + max(clu[[i-1]])
		}
	}

	optfun<-gen.opt.par(M=M,k=k,maxiter=maxiter,approach=approach,trace.iter=trace.iter,save.initial.param = save.initial.param,skip.par=skip.par,save.checked.par=save.checked.par,merge.save.skip.par=merge.save.skip.par,check.skip=check.skip,switch.names=switch.names,...)
	eval(optfun)
	return(opt.par.tmp(M=M,clu=clu,k=k,approach=approach,...))
}

