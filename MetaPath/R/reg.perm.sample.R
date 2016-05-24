reg.perm.sample <-
function(expr,testgroup,nperm) {

## expr=expr; testgroup=testgroup; nperm=nperm;
	obs= cor.func(expr,testgroup)$tt[,1]
	perms= matrix(NA,  nrow(expr),  nperm)
	perms.mtx=matrix(NA,nrow=nperm,ncol=length(testgroup))
	
	for(i in 1:nperm){
		perms.mtx[i,]=sample(1:length(testgroup), size=length(testgroup))
	}

	for(t1 in 1:nperm){
		perms[,t1]=cor.func(expr,testgroup[perms.mtx[t1,]])$tt[,1]
	}
	rownames(perms)=rownames(expr)
	colnames(perms)=paste('B',1:nperm,sep="")
	obs=abs(obs) ### vector
	perms=abs(perms) ## matrix
    return(list(obs=obs, perms=perms,perms.mtx=perms.mtx))
}
