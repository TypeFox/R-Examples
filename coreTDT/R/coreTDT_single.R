#source("/Users/yujiang/Research/compound_heterozygosity/for_R_package_temp/compHet_TDT_single.R")


coreTDT_single = function(ped,maf.threshold=1,qc.proportion=0.8,geneid=NA,control.maf=NULL){
		
	nfamily = dim(ped)[1]/3
	nsnp = dim(ped)[2]
	##variants selection based on maf in evs 
	
	
	
	if(nsnp < 1){
		cat(geneid,nsnp,NA,NA,NA,NA,NA,NA,NA,NA,"\n")
		return(NULL)
	}
	
	
	##check the coverage information
	coverage_infor = apply(ped,2,function(x) mean(apply(matrix(x,ncol=3),1,function(y) any(is.na(y)))))
		
	var_select = which(coverage_infor <= (1-qc.proportion))
	
	#var_select = c(1:nsnp)
	
	nsnp_select = length(var_select)
	
	
	#cat("nsnp_select: ",nsnp_select,"\n")
	
	child.geno.select1 =as.matrix(ped[1:nfamily,var_select])
	parent.geno.select1 = array(NA,dim=c(nfamily,nsnp_select,2))
	parent.geno.select1[1:nfamily,1:nsnp_select,1] = as.matrix(ped[(nfamily+1):(2*nfamily),var_select])
	parent.geno.select1[1:nfamily,1:nsnp_select,2] =as.matrix(ped[(2*nfamily+1):(3*nfamily),var_select])
	
	if(is.null(control.maf)){
		parent_maf = apply(parent.geno.select1,2,function(x) (0.5 * mean(x,na.rm=T)))
	
		var_select2 = which(parent_maf <= maf.threshold)

	}else{
		
		var_select2 = which(control.maf[var_select] <= maf.threshold)
	}	
			
	
	nsnp_select2 = length(var_select2)
	
	#cat("nsnp_select2: ",nsnp_select2,"\n")
	if(nsnp_select2 < 1){
		cat(geneid,nsnp_select2,NA,NA,NA,NA,NA,NA,NA,NA,"\n")
		return(NULL)
	}
	
	child.geno.select2 =array(NA,dim=c(nfamily,nsnp_select2))
	parent.geno.select2 = array(NA,dim=c(nfamily,nsnp_select2,2))
	
	child.geno.select2[1:nfamily,1:nsnp_select2] = child.geno.select1[,var_select2]
	parent.geno.select2[1:nfamily,1:nsnp_select2,1] = parent.geno.select1[,var_select2,1]
	parent.geno.select2[1:nfamily,1:nsnp_select2,2] = parent.geno.select1[,var_select2,2]
	
	
	
	
	
	
	
	##compute the p values
		
	result = compHet_TDT_v6(parent.geno.select2,child.geno.select2)
	
	return(result)
}	






