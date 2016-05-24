`evalRedundancy` <-
function(mat,ia,red=0){
	vec.snps<-strsplit(ia," & ")[[1]]
	n.snps<-length(vec.snps)
	if(n.snps==1)
		return(ia)
	vec.jack<-character(n.snps)
	for(i in 1:n.snps)
		vec.jack[i]<-paste(vec.snps[-i],collapse=" & ")
	vec.red<-numeric(n.snps)
	full<-with(mat, sum(eval(parse(text=ia))))
	for(i in 1:n.snps)
		vec.red[i]<-with(mat, sum(eval(parse(text=vec.jack[i]))))
	if(any(vec.red-full<=red)){
		warning("A redundant SNP in the explanatory interactions is removed.",call.=FALSE)
		ids<-which.min(vec.red-full)[1]
		new.ia<-vec.jack[ids]
		ia<-evalRedundancy(mat,new.ia,red=red)
	}
	ia
}

