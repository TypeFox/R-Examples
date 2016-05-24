"sedist" <-
function(
	M,	#matrix (of a network)
	method="default", 	# the a method used to compute distances - any of the methods alloed by functions dist, cor or cov {all package::stats} or just "cor" or "cov" (given as character)
	fun="default",	#which function should be used to comput distacnes (given as character),
	fun.on.rows="default", # for non-standard function - does it compute measure on rows (such as cor, cov,...) of the data matrix.
#	stats.dist.cor.cov=TRUE,	#call "stats::dist", "stats::cor" or "stats::cov", not "dist", "cor" or "cov", if nonstandard functions are used, they should exemp the same arguments as those in package stats
	handle.interaction="switch",	#how should the interaction between the vertices analysed be handled:
						# "switch" (the default) - assumes that when comparing units i and j, M[i,i] should be compared with M[j,j] and M[i,j] with M[j,i]
						# "switch2" - the same as above, only that each pair occours only once
						# "ignore" (diagonal) - Diagonal is ignored
						# "none" - the matrix is used "as is"
	use = "pairwise.complete.obs",	#for use with methods "cor" and "cov", for other methods (the default option should be used if handle.interaction=="ignore"), "pairwise.complete.obs" are always used, if stats.dist.cor.cov=TRUE
	#p=2	,#The power of the Minkowski distance in functin dist if stats.dist.cor.cov=TRUE
	... #other argumets passed to fun
)
{
	method<-match.arg(method, choices=c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski","pearson", "kendall", "spearman","dist","cor", "cov", "default"))
	if(any(method=="default", fun=="default")){
		if(all(method=="default", fun=="default")){
			fun<-"dist"
			method<-"euclidean"
		} else if(fun=="default"){
			if(method %in% c("pearson", "kendall", "spearman")) fun<-"cor"
			if(method %in% c("cor", "cov")){
				fun<-method
				method<-"pearson"
			}
			if(method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) fun<-"dist"
		} else {
			if(fun %in% c("cor","cov")) method<-"pearson"
			if(fun=="dist") method<-"euclidean"
		}
	}

	if(handle.interaction=="ignore"&& fun %in% c("cor","cov") && !exists(use))warning("The option use='pairwise.complete.obs' should be used with handle.interaction=='ignore' && fun %in% c('cor','cov')")

#	if(fun %in% c("dist", "cor" or "cov") && stats.dist.cor.cov) fun<-paste("stats::",fun,sep="")
	if(fun.on.rows=="default") if(fun %in% c("cor","cov")){
		fun.on.rows<-TRUE
	} else fun.on.rows<-FALSE

	n<-dim(M)[1]
	if(n!=dim(M)[2]) stop("This function is suited for one-mode networks only")
  if(fun %in% c("cor", "cov")) usearg<-list(use=use) else usearg<-NULL #usearg

	if(handle.interaction %in% c("switch","switch2")){
		X<-cbind(M,t(M))
		n<-dim(M)[1]
		res<-matrix(NA,ncol=n,nrow=n)
		for(i in 2:n)for(j in seq(length=(i-1))){
			jind<-seq(length=2*n)
			jind[i]<-j
			jind[j]<-i
			jind[n+i]<-ifelse(handle.interaction=="switch",n+j,NA)
			jind[n+j]<-ifelse(handle.interaction=="switch",n+i,NA)
			Xij<-rbind(X[i,],X[j,jind])
			if(fun.on.rows)Xij<-t(Xij)
			res[i,j]<-do.call(fun,args=c(list(x=Xij, method=method,...),usearg))
		}
		res<-as.dist(res)
	}else{
		if(handle.interaction=="ignore") diag(M)<-NA
		X<-cbind(M,t(M))
		if(fun.on.rows)X<-t(X)
		res<-do.call(fun,args=list(x=X, method=method,...))
	}
	if(class(res)=="dist")attr(res,"Labels")<-rownames(M)
	if(is.matrix(res))dimnames(res)<-dimnames(M)
	return(res)	
}

