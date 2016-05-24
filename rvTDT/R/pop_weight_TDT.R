library(CompQuadForm)

kernel_TDT=function(parent.geno,child.geno,snp.weight){
        
        nfamily = dim(child.geno)[1]
        nsnp = dim(child.geno)[2]
        
        u = matrix(NA,ncol=nsnp,nrow=nfamily)

        ##compute x-E(x|x_p) for each sample at each locus
        for(i in c(1:nsnp)) {
                
                family.hap = cbind(parent.geno[,i,1],parent.geno[,i,2],child.geno[,i])
                u[,i] = apply(family.hap,1,function(x) (x[3]-mean(x[1:2])))
        }
        
        na_pos = which(is.na(u) ==T,arr.ind=T)
		u[na_pos] = 0

        ##compute teh statistics
        
        Q = sum(snp.weight^2 * (apply(u,2,sum))^2) 
        
        ##compute the covariance of WW(X-E(X|X_p))'(Y-mu)
        covariance = matrix(0,ncol=nsnp,nrow=nsnp)
        
        for(i in c(1:nsnp)){
                for(j in c(1:nsnp)) {
                        covariance[i,j] = snp.weight[j]^2 * sum(u[,i] *u[,j])                   
                }               
        }
        
        
	eigens =eigen(covariance,only.values =T)$values
	p = davies(Q,eigens)$Qq

        return(p)
        
        
}



##linear combination test
lc_TDT = function(parent.geno,child.geno,snp.weight) {
	
	
	nfamily = dim(child.geno)[1]
	nsnp = dim(child.geno)[2]
	
	u = matrix(NA,ncol=nsnp,nrow=nfamily)
	

	##compute x-E(x|x_p) for each sample at each locus
	for(i in c(1:nsnp)) {
		family.hap = cbind(parent.geno[,i,1],parent.geno[,i,2],child.geno[,i])
		u[,i] = apply(family.hap,1,function(x) (x[3]-mean(x[1:2])))
	}
	
	na_pos = which(is.na(u) ==T,arr.ind=T)
	u[na_pos] = 0

	
	u.weight = apply(u,1,function(x) sum(x*snp.weight))
	u.sum = sum(u.weight)
	v.emp = sum(u.weight^2)
	s = ifelse(v.emp!=0,u.sum^2/v.emp,0)
	p = 1-pchisq(s,df=1)
	
	return(p)
}


TDT_permutation = function(parent.geno,child.geno,snp.weight1,snp.weight2,snp.weight3,nperm) {
	
	
	nfamily = dim(child.geno)[1]
	nsnp = dim(child.geno)[2]
	
	u = matrix(NA,ncol=nsnp,nrow=nfamily)
	

	##compute x-E(x|x_p) for each sample at each locus
	for(i in c(1:nsnp)) {
		family.hap = cbind(parent.geno[,i,1],parent.geno[,i,2],child.geno[,i])
		u[,i] = apply(family.hap,1,function(x) (x[3]-mean(x[1:2])))
	}
	
	na_pos = which(is.na(u) ==T,arr.ind=T)
	u[na_pos] = 0

	
	u.weight1 = apply(u,1,function(x) sum(x*snp.weight1))
	u.sum1 = sum(u.weight1)	
	Q1 = sum(snp.weight1^2 * (apply(u,2,sum))^2)
	
	u.weight2 = apply(u,1,function(x) sum(x*snp.weight2))
	u.sum2 = sum(u.weight2)	
	Q2 = sum(snp.weight2^2 * (apply(u,2,sum))^2)
	
	u.weight3 = apply(u,1,function(x) sum(x*snp.weight3))
	u.sum3 = sum(u.weight3)		
	Q3 = sum(snp.weight3^2 * (apply(u,2,sum))^2)
	
	u.sum1.perm = rep(NA,nperm)
	u.sum2.perm = rep(NA,nperm)
	u.sum3.perm = rep(NA,nperm)
	
	Q1.perm = rep(NA,nperm)
	Q2.perm = rep(NA,nperm)
	Q3.perm = rep(NA,nperm)
	
	for(iperm in c(1:nperm)){
		set.seed(iperm)
		permPar = sample(c(1,-1),size=nfamily,replace=T)
		u.perm = u*permPar
		u.sum1.perm[iperm] = sum(apply(u.perm,1,function(x) sum(x*snp.weight1)))	
		u.sum2.perm[iperm] = sum(apply(u.perm,1,function(x) sum(x*snp.weight2)))
		u.sum3.perm[iperm] = sum(apply(u.perm,1,function(x) sum(x*snp.weight3)))
		Q1.perm[iperm] = sum(snp.weight1^2 * (apply(u.perm,2,sum))^2)
		Q2.perm[iperm] = sum(snp.weight2^2 * (apply(u.perm,2,sum))^2)
		Q3.perm[iperm] = sum(snp.weight3^2 * (apply(u.perm,2,sum))^2)		
	}
	
	
	p.lc.1 = mean(abs(u.sum1.perm) >= abs(u.sum1))
	p.lc.2 = mean(abs(u.sum2.perm) >= abs(u.sum2))
	p.lc.3 = mean(abs(u.sum3.perm) >= abs(u.sum3))
	p.k.1 = mean(Q1.perm >= Q1)
	p.k.2 =  mean(Q2.perm >= Q2)
	p.k.3 =  mean(Q3.perm >= Q3)
	
	return(c(p.lc.1,p.lc.2,p.lc.3,p.k.1,p.k.2,p.k.3))
}


