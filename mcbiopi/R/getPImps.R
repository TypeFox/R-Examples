getPImps<-function(ltree,type){
	if(type!=1 && (is.na(ltree$coef) || round(ltree$coef,8)==0))
		return(NULL)
	if(type==1 && all(ltree$trees$pick==0))
		return(NULL)
	mat.truth<-generateTruthTab(ltree)
	truth<-ifelse(ltree$coef>=0 | is.na(ltree$coef),1,0)
	ids.truth<-mat.truth[,"outcome"]==truth
	mat.truth<-mat.truth[ids.truth,-ncol(mat.truth),drop=FALSE]
	if(sum(ids.truth)==1)
		vec.primes<-paste(ifelse(mat.truth==1,"","!"),colnames(mat.truth),sep="",
			collapse=" & ")
	else
		vec.primes<-prime.implicants(mat.truth)$vec.primes
	vec.primes
}

