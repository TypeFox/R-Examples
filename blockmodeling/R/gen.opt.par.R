"gen.opt.par" <-
function(
	#function for optimizig partition in blockmodeling
	M,	#matrix
	k,	#number partitions for each mode
#	e1,	#weight of the error of binearized matrix
#	e2,	#weight of the error of valued conenctions
	approach,
	...,	#other arguments to called functions - to 'crit.fun'
	#m,	#suficient value individual cells
	maxiter,	#maximum number of iterations
	mingr=1,	#minimal alowed group size
	maxgr=Inf,	#maximal alowed group size
#	s="default",	#suficient value for colum and row statistics
#	FUN,	#function to calculate row and colum statistics
#	blocks#=c("null","com","reg"),	#permissible block types and their ordering
#	BLOCKS=NULL,	#array of permissible block types and their ordering for all blocks
#	mindim = 2,	#minimal dimension for regulal, dominant and functional blocks
#	save.err.v=FALSE,	#save a vector of errors of all block types for all blocks
	trace.iter,	#save a result of each iteration or only the best (minimal error)
	switch.names,	#should partitions that only differ in group names be considert equal (is c(1,1,2)==c(2,2,1))
	save.initial.param,	#should the initial parameters be saved
	skip.par,	#the partions that are not allowed or were already checked and should be skiped
	save.checked.par,	#should the checked partitions be saved
	merge.save.skip.par, #should the checked partitions be merged with skiped ones
	parOK = NULL, # a function that takes the partition as the argument and retutns TRUE, it partition is ok (satisfies the all conditions)
	parOKaddParam = NULL, # additional parameters for function "partOK"
	check.skip,	#when should the check be preformed:
								# "all"  - before every call to 'crit.fun'
								# "iter" - at the end of eack iteratiton
								# "opt.par"  - before every call to 'opt.par', implemented in opt.these.par and opt.random.par
								# "never" - never
	print.iter=TRUE,
	check.switch=TRUE,
	check.all=TRUE,
	use.for.opt=TRUE
#	force.fun=NULL, #select the function used to evaluate partition
){
	dots<-list(...)
	nmode<-length(k)
	notUseNormMto2rel<-is.null(dots$normMto2rel)||!dots$normMto2rel
	if(approach=="ss"&&all(dots$blocks=="com")&& all(dots$BLOCKS=="com")&&nmode==1&&use.for.opt&&notUseNormMto2rel){
  	fun<-c("opt.par.tmp<-function(M=M,clu=clu,diag=TRUE,...){
			maxiter<-",maxiter,"
			n<-dim(M)[1]
			clu<-as.integer(factor(clu))
			k<-max(clu)
			M[,]<-as.double(M)
			res<-.Fortran('optparsscom',M=M,clu=clu,diag=diag,maxiter=as.integer(maxiter),n=as.integer(dim(M)[1]),k=k,err=as.double(0),E=diag(k)*as.double(0),BM=diag(k)*as.double(0),cluM=matrix(as.integer(0),nrow=50,ncol=dim(M)[1]),nbest=as.integer(0),iter=as.integer(0),printIter=",print.iter,")
			best1<-res[c(2,7,8,9)]
      IM<-best1$E
      IM[,] <- 'com'
      best1$IM<-IM
			res<-list(M=res$M,best=c(list(best1=best1), if(res$nbest>1) apply(res$cluM[seq(res$nbest),],1,function(x)list(clu=x))[-1]else NULL),nIter=res$iter)
			class(res)<-'opt.par'
			return(res)
		}")
#		cat(fun,file="tmp.R")
		return(parse(text=fun))
	}

	if(approach=="ss"&&all(dots$blocks=="com")&& all(dots$BLOCKS=="com")&&nmode==2&&use.for.opt&&notUseNormMto2rel&&length(dim(M))==2){
  	fun<-c("opt.par.tmp<-function(M=M,clu=clu,diag=TRUE,...){
			maxiter<-",maxiter,"
			n1<-dim(M)[1]
			n2<-dim(M)[2]
			clu1<-as.integer(factor(clu[[1]]))
			clu2<-as.integer(factor(clu[[2]]))
			k1<-max(clu1)
			k2<-max(clu2)
			M[,]<-as.double(M)
			res<-.Fortran('optparsscomtm',M=M,clu1=clu1,clu2=clu2,maxiter=as.integer(maxiter),n1=as.integer(dim(M)[1]),n2=as.integer(dim(M)[2]),k1=k1,k2=k2,err=as.double(0),E=matrix(as.double(0),nrow=k1,ncol=k2),BM=matrix(as.double(0),nrow=k1,ncol=k2),cluM1=matrix(as.integer(0),nrow=50,ncol=dim(M)[1]),cluM2=matrix(as.integer(0),nrow=50,ncol=dim(M)[2]),nbest=as.integer(0),iter=as.integer(0),printIter=",print.iter,")
			best1<-res[c('err','E','BM')]
			best1$clu<-list(res$clu1,res$clu2)
			IM<-best1$E
			IM[,] <- 'com'
			best1$IM<-IM
			bestClus<-list()
			if(res$nbest>1)for(i in 1:res$nbest){
				bestClus<-c(bestClus,list(list(clu=list(res$cluM1[i,],res$cluM2[i,]))))
			}
			res<-list(M=res$M,best=c(list(best1=best1), if(res$nbest>1) bestClus),nIter=res$iter)
			class(res)<-'opt.par'
			return(res)
		}")
		fun<-paste(fun,collapse="")
#		cat(fun,file="tmp.R")
		return(parse(text=fun))
	}

	if(approach=="ss"&&all(dots$blocks=="com")&& all(dots$BLOCKS=="com")&&nmode==2&&use.for.opt&&notUseNormMto2rel&&length(dim(M))==3){
  	fun<-c("opt.par.tmp<-function(M=M,clu=clu,diag=TRUE,...){
			maxiter<-",maxiter,"
			n1<-dim(M)[1]
			n2<-dim(M)[2]
			nr<-dim(M)[3]
			clu1<-as.integer(factor(clu[[1]]))
			clu2<-as.integer(factor(clu[[2]]))
			k1<-max(clu1)
			k2<-max(clu2)
			M[,,]<-as.double(M)
			res<-.Fortran('optparsscomtmmorerel',M=M,clu1=clu1,clu2=clu2,maxiter=as.integer(maxiter),nr=as.integer(dim(M)[3]),n1=as.integer(dim(M)[1]),n2=as.integer(dim(M)[2]),k1=k1,k2=k2,err=as.double(0),E=matrix(as.double(0),nrow=k1,ncol=k2),BM=array(as.double(0),dim=c(k1,k2,nr)),cluM1=matrix(as.integer(0),nrow=50,ncol=dim(M)[1]),cluM2=matrix(as.integer(0),nrow=50,ncol=dim(M)[2]),nbest=as.integer(0),iter=as.integer(0),printIter=",print.iter,")
			best1<-res[c('err','E','BM')]
			best1$clu<-list(res$clu1,res$clu2)
			IM<-best1$E
			IM[,] <- 'com'
			best1$IM<-IM
			bestClus<-list()
			if(res$nbest>1)for(i in 1:res$nbest){
				bestClus<-c(bestClus,list(list(clu=list(res$cluM1[i,],res$cluM2[i,]))))
			}
			res<-list(M=res$M,best=c(list(best1=best1), if(res$nbest>1) bestClus),nIter=res$iter)
			class(res)<-'opt.par'
			return(res)
		}")
		fun<-paste(fun,collapse="")
#		cat(fun,file="tmp.R")
		return(parse(text=fun))
	}




	fun<-c("opt.par.tmp<-function(M=M,clu=clu,k=k,approach=approach,",if(!is.null(parOK)) "parOK=parOK, parOKaddParam = parOKaddParam,","...){\n")
	if(save.initial.param) fun<-c(fun,"initial.param<-tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return('error'))\n")	#saves the inital parameters

	fun<-c(fun,"	critfun<-gen.crit.fun(M=M,k=k,approach=approach,changeT=TRUE,...)

	eval(critfun$fun1)\n")

	if(trace.iter) {
		fun<-c(fun,"best.crit.iter<-list(NULL)\n")
	}#else best.crit.iter<-NULL


	fun<-c(fun,"best.crit<-list(NULL)
	
	res.last.iter<-best.crit[[1]]<-crit.fun.tmp(M=useM,clu)
	eval(critfun$fun2)
", if(save.checked.par) "checked.par<-list(clu)\n" else NULL,
"	iter<-0
	imp<-TRUE

	cat('Starting partition:',",if(nmode>1) "unlist(" else NULL, "clu",if(nmode>1) ")" else NULL, ",'\\n')
	cat('Starting error:',best.crit[[1]]$err,'\\n\\n')


	while((iter<maxiter)&imp){
		imp<-FALSE
		iter<-iter+1
", "		", if(nmode==1) "tclu<-table(clu)" else "tclu<-lapply(clu,table)","\n",
if(nmode>1)c("	for(imode in 1:nmode)",if(!check.all)"if(!imp)" else NULL,"{\n") else NULL,
"		for(i1 in ",if(nmode>2)"cumsum(c(0,k)[1:imode])+" else NULL,"1:k",if(nmode>1) "[[imode]]" else NULL, ")",if(!check.all)"if(!imp)" else NULL,"{
			for(i2 in setdiff(",if(nmode>2)"cumsum(c(0,k)[1:imode])+" else NULL,"1:k",if(nmode>1) "[[imode]]" else NULL,",i1))",if(!check.all)"if(!imp)" else NULL,"{
				for(i3 in 1:tclu",if(nmode>1) "[[imode]]" else NULL,"[i1])",if(!check.all)"if(!imp)" else NULL,"{
					if(tclu",if(nmode>1) "[[imode]]" else NULL,"[i1]>",mingr,"&tclu",if(nmode>1) "[[imode]]" else NULL,"[i2]<",maxgr,") {	#checkin if movement to another group improves criteria function
						tempclu<-clu
						tempclu",if(nmode>1) "[[imode]]" else NULL, "[tempclu",if(nmode>1) "[[imode]]" else NULL, "==i1][i3]<-i2
", if(!is.null(parOK))"if(parOK(tempclu,parOKaddParam)){" else NULL ,if(check.skip=="all") c("							if(
								ifelse(
									is.null(skip.par),
									TRUE,",
									if(switch.names) c("
									!any(sapply(skip.par,",if(nmode>1) "function(x,clu2)rand2(unlist(x),clu2)" else "rand2", ",clu2=",if(nmode>1) "unlist(" else NULL, "tempclu",if(nmode>1) ")" else NULL, ")==1)"
									)else c("
									!any(sapply(skip.par,function(x)all(",if(nmode>1) "unlist(" else NULL, "x",if(nmode>1) ")" else NULL, "==",if(nmode>1) "unlist(" else NULL, "tempclu",if(nmode>1) ")" else NULL, ")))"
									),"
								)
							) #checking if the partition should be ignored
							{
") else NULL,"						temp.crit<-crit.fun.tmp(M=useM,clu=tempclu,change=c(i1,i2)",if(nmode==2) ", modechange=imode" else NULL, ",res.old=res.last.iter)

						if(temp.crit$err<best.crit[[1]]$err) {
							best.crit<-list(NULL)
							best.crit[[1]]<-temp.crit
							imp<-TRUE
						} else if(temp.crit$err==best.crit[[1]]$err){
							unique<-TRUE
							for(i in 1:length(best.crit)){
								",if(switch.names){c(
								"if(rand(table(",if(nmode>1) "unlist(" else NULL, "best.crit[[i]]$clu",if(nmode>1) ")" else NULL, ",",if(nmode>1) "unlist(" else NULL, "temp.crit$clu",if(nmode>1) ")" else NULL, "))==1) unique<-FALSE")
								}else {c("
								if(all(",if(nmode>1) "unlist(" else NULL, "best.crit[[i]]$clu",if(nmode>1) ")" else NULL, "==",if(nmode>1) "unlist(" else NULL, "temp.crit$clu",if(nmode>1) ")" else NULL, ")) unique<-FALSE")
								},"
							}
							if(unique) {best.crit[[length(best.crit)+1]]<-temp.crit}
						}
", if(check.skip=="all") "							}
" else NULL,if(!is.null(parOK))"}" else NULL,"					}
					",if(check.switch)c("if(i1<i3){	#checkin if exchange of two units (from different groups) improves criteria function
						for(i4 in 1:tclu",if(nmode>1) "[[imode]]" else NULL,"[i2])",if(!check.all)"if(!imp)" else NULL,"{
							tempclu<-clu
							tempclu",if(nmode>1) "[[imode]]" else NULL,"[clu",if(nmode>1) "[[imode]]" else NULL,"==i1][i3]<-i2
							tempclu",if(nmode>1) "[[imode]]" else NULL,"[clu",if(nmode>1) "[[imode]]" else NULL,"==i2][i4]<-i1
", if(!is.null(parOK))"if(parOK(tempclu,parOKaddParam)){" else NULL , if(check.skip=="all") c("							if(
								ifelse(
									is.null(skip.par),
									TRUE,",
									if(switch.names)c("
									!any(sapply(skip.par,",if(nmode>1) "function(x,clu2)rand2(unlist(x),clu2)" else "rand2", ",clu2=",if(nmode>1) "unlist(" else NULL, "tempclu",if(nmode>1) ")" else NULL, ")==1)"
									)else c("
									!any(sapply(skip.par,function(x)all(",if(nmode>1) "unlist(" else NULL, "x",if(nmode>1) ")" else NULL, "==",if(nmode>1) "unlist(" else NULL, "tempclu",if(nmode>1) ")" else NULL, ")))"
									),"
								)
							) #checking if the partition should be ignored
							{
") else NULL,"							temp.crit<-crit.fun.tmp(M=useM,clu=tempclu,change=c(i1,i2)",if(nmode==2) ", modechange=imode" else NULL, ",res.old=res.last.iter)
							if(temp.crit$err<best.crit[[1]]$err) {
								best.crit<-list(NULL)
								best.crit[[1]]<-temp.crit
								imp<-TRUE
							} else if(temp.crit$err==best.crit[[1]]$err){
								unique<-TRUE
								for(i in 1:length(best.crit)){"
									,if(switch.names){c("
									if(rand(table(",if(nmode>1) "unlist(" else NULL, "best.crit[[i]]$clu",if(nmode>1) ")" else NULL, ",",if(nmode>1) "unlist(" else NULL, "temp.crit$clu",if(nmode>1) ")" else NULL, "))==1) unique<-FALSE")
									}else {c("
									if(all(",if(nmode>1) "unlist(" else NULL, "best.crit[[i]]$clu",if(nmode>1) ")" else NULL, "==",if(nmode>1) "unlist(" else NULL, "temp.crit$clu",if(nmode>1) ")" else NULL, ")) unique<-FALSE")
									},"
								}
								if(unique) {best.crit[[length(best.crit)+1]]<-temp.crit}
							}
", if(check.skip=="all") "							}
" else NULL,if(!is.null(parOK))"}" else NULL,"						}
					}") else NULL,"
				}
			}
		}
	",if(nmode>1) "}\n\t" else NULL, if(trace.iter) "best.crit.iter[[iter]]<-best.crit" else NULL,"
	res.last.iter<-best.crit[[1]]
	clu<-best.crit[[1]]$clu
", if(check.skip=="iter") c("		if(
			ifelse(
				is.null(skip.par),
				FALSE,",
				if(switch.names){c("
				any(sapply(skip.par,",if(nmode>1) "function(x,clu2)rand2(unlist(x),clu2)" else "rand2", ",clu2=",if(nmode>1) "unlist(" else NULL, "clu",if(nmode>1) ")" else NULL, ")==1)")
				}else{c("
				any(sapply(skip.par,function(x)all(",if(nmode>1) "unlist(" else NULL, "x",if(nmode>1) ")" else NULL, "==",if(nmode>1) "unlist(" else NULL, "clu",if(nmode>1) ")" else NULL, ")))")
				},"
			)
		) #checking if the partition should be ignored
		{
			imp<-FALSE
			best.crit[[1]]$err<-NA
			cat('Optimization ended - partition was in \\'skip.par\\'\\n')
		}
") else NULL, if(save.checked.par) "	if(imp) checked.par<-c(checked.par,list(clu))
" else NULL, if(print.iter) c("	cat('End of iteration',iter,'\\n')
	cat('Current partition: ',",if(nmode>1) "unlist(" else NULL, "clu",if(nmode>1) ")" else NULL, ",'\\n')
	cat('Current error:',best.crit[[1]]$err,'\\n\\n')	#the end of an iteration
") else NULL,"	}
", if(merge.save.skip.par && save.checked.par) "	checked.par<-c(skip.par, checked.par)
" else NULL,"	cat('Function completed\\n')
	cat('Final (1st) paritition:',",if(nmode>1) "unlist(" else NULL,"best.crit[[1]]$clu",if(nmode>1) ")" else NULL,",'\\n')
	cat(length(best.crit),'solution(s) with minimal error =', best.crit[[1]]$err, 'found.','\\n')

	res<-list(M=M,best=best.crit", if(trace.iter) ",iter=best.crit.iter" else NULL,",call=match.call()",if(save.initial.param) ",initial.param = initial.param" else NULL, if(save.checked.par) ", checked.par=checked.par" else NULL,",nIter=iter)
	class(res)<-'opt.par",if(nmode>1) ".mode" else NULL, "'
	return(res)
}\n")
#cat(fun,sep="",file="tmp.R")
return(parse(text=paste(fun,sep="",collapse="")))
}

