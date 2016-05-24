
rvTDT = function(ped,evs,maf.threshold = 1, qc.proportion =0.8){


#	pedfile = paste(geneid,".ped2",sep="")
#	mapfile = paste(geneid,".map2",sep="")
#	evsfile = paste(geneid,"_evs.map2",sep="")

#	ped = try(read.table(pedfile),silent=T)
#	map = try(read.table(mapfile),silent=T)
#	evs = try(read.table(evsfile,row.name=1),silent=T)


	nfamily = dim(ped)[1]/3
	nsnp = dim(ped)[2]
	##variants selection based on maf in evs 
	
	result=list()
	class(result)="rvTDT"

#	result$geneid = geneid
	result$nsnptot = nsnp
	result$nfamily = nfamily
	
	if(nsnp <= 1){
		result$nsnpcompute = nsnp
		result$p_lc_1 = NA
		result$p_lc_maf = NA
		result$p_lc_pc = NA
		result$p_k_1 = NA
		result$p_k_maf = NA
		result$p_k_pc = NA
		return(result)
	}

	##check the NA in files
	
	coverage_infor = apply(ped,2,function(x) mean(apply(matrix(x,ncol=3),1,function(y) any(is.na(y)))))
	
	var_select = which(coverage_infor <= (1-qc.proportion))
	
#	var_select = c(1:nsnp)

	nsnp_select = length(var_select)
	
	if(nsnp_select <=1){
		result$nsnpcompute = nsnp_select
		result$p_lc_1 = NA
		result$p_lc_maf = NA
		result$p_lc_pc = NA
		result$p_k_1 = NA
		result$p_k_maf = NA
		result$p_k_pc = NA
		return(result)

	}
	
	

	child.geno.select1 =as.matrix(ped[1:nfamily,var_select])
	parent.geno.select1 = array(NA,dim=c(nfamily,nsnp_select,2))
	parent.geno.select1[1:nfamily,1:nsnp_select,1] = as.matrix(ped[(nfamily+1):(2*nfamily),var_select])
	parent.geno.select1[1:nfamily,1:nsnp_select,2] =as.matrix(ped[(2*nfamily+1):(3*nfamily),var_select])

	evs_sum = as.matrix(evs[var_select,1:3])
	parent_sum = matrix(NA,nrow=nsnp_select,ncol=3)

	parent_sum[,1] = apply(parent.geno.select1,2,function(x) sum(x==2,na.rm=T))
	parent_sum[,2] = apply(parent.geno.select1,2,function(x) sum(x==1,na.rm=T))
	parent_sum[,3] = apply(parent.geno.select1,2,function(x) sum(x==0,na.rm=T))

	parent_maf = apply(parent.geno.select1,2,function(x) (0.5 * mean(x,na.rm=T)))
##compute weights

	unweighted = rep(1,nsnp_select)
	pc_weights = TrendWeights(evs_sum,parent_sum)
	maf_weights = dbeta(parent_maf,1,25)


##evs_maf, set weight for common variants to be 0
	evs_maf = apply(evs_sum,1, function(x) ((x[1]*2 + x[2])/(sum(x)*2)))
#	unweighted[which(evs_maf > maf.threshold)] = 0
	pc_weights[which(evs_maf > maf.threshold)] = 0
#	maf_weights[which(evs_maf > maf.threshold)] = 0

	unweighted[which(parent_maf > maf.threshold)] = 0
	maf_weights[which(parent_maf > maf.threshold)] = 0

	#unweighted[which(parent_maf ==0)] = 0
	#pc_weights[which(parent_maf ==0)] = 0
	#maf_weights[which(parent_maf ==0)] = 0

	#unweighted[which(parent_maf > 0.05)] = 0
	#pc_weights[which(parent_maf > 0.05)] = 0
	#maf_weights[which(parent_maf > 0.05)] = 0

	result$nsnpcompute = length(which(unweighted !=0))


##compute the p values

	T_lc_1  = lc_TDT(parent.geno.select1,child.geno.select1,unweighted)
	T_k_1 = kernel_TDT(parent.geno.select1,child.geno.select1,unweighted)

	T_lc_maf  = lc_TDT(parent.geno.select1,child.geno.select1,maf_weights)
	T_k_maf = kernel_TDT(parent.geno.select1,child.geno.select1,maf_weights)

	T_lc_pc  = lc_TDT(parent.geno.select1,child.geno.select1,pc_weights)
	T_k_pc = kernel_TDT(parent.geno.select1,child.geno.select1,pc_weights)

	result$p_lc_1 = T_lc_1
	result$p_lc_maf = T_lc_maf
	result$p_lc_pc = T_lc_pc
	result$p_k_1 = T_k_1
	result$p_k_maf = T_k_maf
	result$p_k_pc = T_k_pc

	
	return(result)


}



