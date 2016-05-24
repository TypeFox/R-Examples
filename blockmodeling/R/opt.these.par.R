"opt.these.par" <-
function(
	M,	#matrix (network)
	partitions,	#a list of partitions to check
	approach,
	...,
	return.all=FALSE,	#if 'FALSE', solution for only the best (one or more) partition/s is/are returned
	return.err=TRUE,	#Should a vector of errors be returned
	skip.allready.checked.par=TRUE,	#if 'TRUE',the partitions that were already checked when runing 'opt.par' form different statrting points will be skiped
	maxiter=50,	#maximum number of iterations
	#m=NULL,	#suficient value individual cells
	#cut=min(M[M>0]),   #
	#BLOCKS=NULL,	#array of permissible block types and their ordering for all blocks
	trace.iter=FALSE,	#save a result of each iteration or only the best (minimal error) (an argument of "gen.opt.par")
	switch.names=NULL,	#should partitions that only differ in group names be considert equal (is c(1,1,2)==c(2,2,1))
	save.initial.param=TRUE,	#should the initial parameters be saved
	skip.par=NULL,	#the partions that are not allowed or were already checked and should be skiped
	save.checked.par=!is.null(skip.par),	#should the checked partitions be saved
	merge.save.skip.par=all(!is.null(skip.par),save.checked.par), #should the checked partitions be merged with skiped ones
	check.skip="never",	#when should the check be preformed:
								# "all"  - before every call to 'crit.fun'
								# "iter" - at the end of eack iteratiton
								# "opt.par"  - before every call to 'opt.par', implemented in opt.these.par and opt.random.par
								# "never" - never
	print.iter=FALSE

){
	if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
	
	res<-list(NULL)
	err<-NULL
	nIter<-NULL

	dots<-list(...)
	if(is.null(switch.names)){
		switch.names<-is.null(dots$BLOCKS)
	}
	
	clu<-partitions[[1]]
	if(is.list(clu)){
		k<-sapply(clu,function(x)length(unique(x)))
	} else k<-length(unique(clu))

	nmode<-length(k)
	
  	optfun<-gen.opt.par(M=M,k=k,maxiter=maxiter,approach=approach,switch.names=switch.names,trace.iter=trace.iter,save.initial.param = save.initial.param,skip.par=skip.par,save.checked.par=save.checked.par,merge.save.skip.par=merge.save.skip.par,check.skip=check.skip,print.iter=print.iter,...)

	on.exit({
		res1 <- res[which(err==min(err, na.rm = TRUE))]
	
		best<-NULL
		best.clu<-NULL
		for(i in 1:length(res1)){
			for(j in 1:length(res1[[i]]$best)){
				if(
					ifelse(is.null(best.clu),
						TRUE,
						if(nmode==1) ifelse(switch.names,
							!any(sapply(best.clu,rand2,clu2=res1[[i]]$best[[j]]$clu)==1),
							!any(sapply(best.clu,function(x)all(x==res1[[i]]$best[[j]]$clu)))
						) else ifelse(switch.names,
							!any(sapply(best.clu,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(res1[[i]]$best[[j]]$clu))==1),
							!any(sapply(best.clu,function(x)all(unlist(x)==unlist(res1[[i]]$best[[j]]$clu))))
						)
					)		
				){ 
					best<-c(best,res1[[i]]$best[j])
					best.clu<-c(best.clu,list(res1[[i]]$best[[j]]$clu))
				}
			}
		}
	
		names(best)<-paste("best",1:length(best),sep="")
	
		cat("\n\nOptimization of all partitions completed\n")
		cat(length(best),"solution(s) with minimal error =", min(err,na.rm=TRUE), "found.","\n")
	
	
		checked.par<-list(checked.par=skip.par)
		call<-list(call=match.call())
		best<-list(best=best)
		if(return.all) res<-list(res=res) else res<-NULL
		if(return.err) err<-list(err=err) else err<-NULL
		if(!exists("initial.param")){
			initial.param<-NULL
		} else initial.param=list(initial.param)
		
		res<-c(list(M=M),res,best,err,list(nIter=nIter),checked.par,call,initial.param=initial.param)
		class(res)<-"opt.more.par"
		return(res)
	})


	eval(optfun)
	npar<-length(partitions)
	for(i in 1:npar){
		cat("\n\nStarting optimization of the partiton",i,"of",npar,"partitions.\n")
		if(
			ifelse(check.skip=="never",
				TRUE,
				ifelse(is.null(skip.par),
					TRUE,
					if(nmode==1) ifelse(switch.names,
						!any(sapply(skip.par,rand2,clu2=res1[[i]]$best[[j]]$clu)==1),
						!any(sapply(skip.par,function(x)all(x==res1[[i]]$best[[j]]$clu)))
					) else ifelse(switch.names,
						!any(sapply(skip.par,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(res1[[i]]$best[[j]]$clu))==1),
						!any(sapply(skip.par,function(x)all(unlist(x)==unlist(res1[[i]]$best[[j]]$clu))))
					)
				)
			)
		) #checking if the partition should be ignored
		{	

			res[[i]]<-opt.par.tmp(
					M=M,
					clu=partitions[[i]],
					k=k,
					approach=approach,
					m=m,
					skip.par=skip.par,
					...
			)
			
			err[i]<-res[[i]]$best[[1]]$err
			nIter[i]<-res[[i]]$nIter
			if(skip.allready.checked.par) skip.par<-c(skip.par,res[[i]]$checked.par)
		} else cat("Optimization of the partition skipped\n")
	}
}

