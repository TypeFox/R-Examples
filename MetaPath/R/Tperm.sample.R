Tperm.sample <-
function(x, fac, nperm) {
    obs = genefilter::rowttests(x, fac, tstatOnly= T)$statistic
	names(obs)=rownames(x)
	perms= matrix(NA,  nrow(x),  nperm)
	
	for(t1 in 1:nperm){
		perms[,t1]=genefilter::rowttests(x, sample(fac), tstatOnly= T)$statistic
	}
	rownames(perms)=rownames(x)
	colnames(perms)=paste('B',1:nperm,sep="")
	
	obs=abs(obs); perms=abs(perms)
		
    return(list(obs=obs, perms=perms))
}
