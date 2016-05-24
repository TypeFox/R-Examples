F.perm.sample <-
function(x, fac, nperm) {
    obs = genefilter::rowFtests(x, fac, var.equal = TRUE)$statistic
	names(obs)=rownames(x)
	perms= matrix(NA,  nrow(x),  nperm)
	
	for(t1 in 1:nperm){
		perms[,t1]=genefilter::rowFtests(x, sample(fac), var.equal = TRUE)$statistic
	}
	rownames(perms)=rownames(x)
	colnames(perms)=paste('B',1:nperm,sep="")
	
	obs=abs(obs); perms=abs(perms)
		
    return(list(obs=obs, perms=perms))
}
