##### internal functions from Dr. Tibshirani's software package GSA
cox.perm.sample <-
function(expr,testgroup,censoring,nperm) {

##  expr=x;testgroup=testlabel;censoring=censoring;nperm=nperm
	obs= coxfunc(expr,testgroup,censoring)$tt
	perms= matrix(NA,  nrow(expr),  nperm)
	perms.mtx=matrix(NA,nrow=nperm,ncol=length(testgroup))
	
	for(i in 1:nperm){
		perms.mtx[i,]=sample(1:length(testgroup), size=length(testgroup))
	}

	for(t1 in 1:nperm){
		perms[,t1]=coxfunc(expr,testgroup[perms.mtx[t1,]],censoring)$tt
	}
	rownames(perms)=rownames(expr)
	colnames(perms)=paste('B',1:nperm,sep="")
	obs=abs(obs) ### vector
	perms=abs(perms) ## matrix
    return(list(obs=obs, perms=perms,perms.mtx=perms.mtx))
}
