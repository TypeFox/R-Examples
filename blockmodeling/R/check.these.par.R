"check.these.par" <-
function(	#saves the resoult of call a to crit.fun of the best partition and only errors for the rest
	M,	#matrix (network)
	partitions,	#partitions to check
	approach, #
	return.err=TRUE,	#if 'FALSE', only the resoults of crit.fun are returned (a list of all (best) soulutions including errors), else the resoult is list
	save.initial.param=TRUE,	#should the initial parameters be saved
#	use.for=TRUE, #should fortran rutines be used when possible
 	force.fun=NULL, #select the function used to evaluate partition
 	...	#paremeters for "gen.crit.fun"
){
	if(save.initial.param){
		initial.param<-list(initial.param=tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")))	#saves the inital parameters
	}else initial.param<-NULL

	err<-NULL
	
	if(is.list(partitions[[1]])){
		k<-sapply(partitions[[1]],function(x)length(unique(x)))
	} else k<-length(unique(partitions[[1]]))
	eval(gen.crit.fun(M=M,k=k,approach=approach,changeT=FALSE,...))


	res<-crit.fun.tmp(M=useM,clu=partitions[[1]])
	best<-list(res)
	err[1]<-res$err
	min.err<-res$err
	 
	for(i in 2:length(partitions)){
		res<-crit.fun.tmp(M=useM,clu=partitions[[i]])
		if(res$err<=min.err){
			if(res$err<min.err) {
				best<-list(res)
				min.err<-res$err
			}else best<-c(best,list(res))
		}
		err[i]<-res$err
	}
	
	names(best)<-paste("best",1:length(best),sep="")
	call<-list(call=match.call())
	best<-list(best=best)
	if(return.err) err<-list(err=err) else err<-NULL
	res<-c(list(M=M),best,err,call,initial.param)
	class(res)<-"check.these.par"
	return(res)
}

