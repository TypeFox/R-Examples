"crit.fun" <-
function(
	#function for generting a function for computing criteria function of a blockmodel and preparing data
	M,	#matrix,
	clu,	#partition
#	e1="default",	#weight of the error of binearized matrix
#	e2="default",	#weight of the error of valued conenctions
	approach,	#the approach used - can be one of "ss","ad","bin","val", "bv", "imp", "bi"
#	cut = min(M[M>0]),
#	m=ifelse(e2==0,1,"default"),	#suficient value for individual cells
#	s="default",	#suficient value for colum and row statistics
#	FUN="max",	#function to calculate row and colum statistics
#	norm=FALSE,
#	blocks=c("null","com","reg"),	#permissible block types and their ordering, can be also on of 'structural', 'regular', 'regular.ext' or 'all'
#	BLOCKS=NULL,	#prespecified model
#	mindim = 2,	#minimal dimension for regulal, dominant and functional blocks
#	block.weights=c(null=1,com=1,rdo=1,cdo=1,reg=1,rre=1,cre=1),	#weights for all block types - if only some of them are specified, the other remain 0
#	fixed.clu=NULL,	#used to specify groups that are not to be evaluated - the blocks that are dependent only on these groups will not be evaluetad
#	initial.IM=NULL,	#an image matrix, usualy the resoult of previou crit.fun call, component "IM" - used only if part of the matrix is fixed - for the solution of the fixed part
#	initial.E=NULL,	#an error matrix, usualy the resoult of previou crit.fun call, component "E" - used only if part of the matrix is fixed - for the solution of the fixed part
#	save.err.v=FALSE,	#save a vector of errors of all block tipes for all blocks
#	max.con.val="non",	#should the largest values be cencored, limited to (larger values set to) - resonoble values are:
#																# "m" or "s" - the maximum value equals the value of the parameter m/s
#																# "non" - no transformation is done
#																# other values larger then parameters m and s and lower the the maximum
#	dn=mean(M>=cut),
  ...
){
	if(is.list(clu)){
		k<-sapply(clu,function(x)length(unique(x)))
		clu<-lapply(clu,function(x)as.integer(factor(x)))
		if(length(k)>2) {
			for(i in 2:length(clu)){
				clu[[i]]<-clu[[i]] + max(clu[[i-1]])
	  		}
	  	}
	} else {
		k<-length(unique(clu))
		clu<-as.integer(factor(clu))
	}
	
	eval(gen.crit.fun(M=M,k=k,approach=approach,changeT=FALSE,...))
	res<-c(list(M=M),crit.fun.tmp(M=useM,clu=clu),call=match.call())
	class(res)<-"crit.fun"
	return(res)
}

