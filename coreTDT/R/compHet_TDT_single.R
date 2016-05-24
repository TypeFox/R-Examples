##features:
##1. allow mendelian error & missing data
##2. remove confusing samples

#source("/Users/yujiang/Research/compound_heterozygosity/for_R_package_temp/phasing_codes_3_25_2014.R")
#source("/Users/yujiang/Research/compound_heterozygosity/for_R_package_temp/pvalue_calculator.R")


compHet_TDT_v6 = function(parent.geno,child.geno){
	
	nfamily = dim(child.geno)[1]
	nsnp = dim(child.geno)[2]

	##recode genotype 
	parent.geno.recode = matrix(NA,nrow=nfamily,ncol=2)
	child.geno.recode = rep(NA,nfamily)
	
	N12 = 0
	N11 = 0
	N112 = 0
	N122 = 0
	nMedErr = 0
	nmissing = 0
	for(ifamily in c(1:nfamily)){
	
		##screen for mendelian error and NA (low coverage families)
		rm.loci = rep(0,nsnp)
		
		for(isnv in c(1:nsnp)){
			##screen for NAs
			if(any(is.na(c(parent.geno[ifamily,isnv,],child.geno[ifamily,isnv])))){
				rm.loci[isnv] = 1
				nmissing = nmissing +1
				next
			}
			if( (any(parent.geno[ifamily,isnv,] ==2) & child.geno[ifamily,isnv] < 1) | (any(parent.geno[ifamily,isnv,] ==0) & child.geno[ifamily,isnv] > 1) | (all(parent.geno[ifamily,isnv,] ==2) & child.geno[ifamily,isnv] < 2) |(all(parent.geno[ifamily,isnv,] ==0) & child.geno[ifamily,isnv] > 0) ){
				rm.loci[isnv] = 1
				nMedErr = nMedErr +1
#				cat(ifamily,isnv,"\n")
			}
			
		}
		
		rm.loci.pos = which(rm.loci ==1)
		
		if(length(rm.loci.pos) == nsnp){
			parent.geno.recode[ifamily,1] = 0
			parent.geno.recode[ifamily,2] = 0
			child.geno.recode[ifamily] = 0
		
		}else if((nsnp - length(rm.loci.pos)) == 1){
			if(nsnp ==1){
				tmp = FamilyPhaseIII(matrix(parent.geno[ifamily,,],ncol=2),child.geno[ifamily,])
				parent.geno.recode[ifamily,1] = tmp[1]
				parent.geno.recode[ifamily,2] = tmp[2]
				child.geno.recode[ifamily] = tmp[3]
				
			}else{
				tmp = FamilyPhaseIII(matrix(parent.geno[ifamily,-rm.loci.pos,],ncol=2),child.geno[ifamily,-rm.loci.pos])
				parent.geno.recode[ifamily,1] = tmp[1]
				parent.geno.recode[ifamily,2] = tmp[2]
				child.geno.recode[ifamily] = tmp[3]
			}

			
		}else if(length(rm.loci.pos) == 0 & nsnp > 1){
			tmp = FamilyPhaseIII(as.matrix(parent.geno[ifamily,,]),child.geno[ifamily,])
			parent.geno.recode[ifamily,1] = tmp[1]
			parent.geno.recode[ifamily,2] = tmp[2]
			child.geno.recode[ifamily] = tmp[3]

		}else{
			tmp = FamilyPhaseIII(as.matrix(parent.geno[ifamily,-rm.loci.pos,]),child.geno[ifamily,-rm.loci.pos])
			parent.geno.recode[ifamily,1] = tmp[1]
			parent.geno.recode[ifamily,2] = tmp[2]
			child.geno.recode[ifamily] = tmp[3]
		}
		
		
		
		if(1 %in% parent.geno.recode[ifamily,]){
			if(all(parent.geno.recode[ifamily,] ==1)){
				N11 =N11+1
				if(child.geno.recode[ifamily] == 2){
					N112 = N112 +1
				}
				
			}else if(2 %in% parent.geno.recode[ifamily,]){
				N12 = N12 +1
				if(child.geno.recode[ifamily] == 2){
					N122 = N122 +1
				}
			}			
		}
	
	}
	
	##compute the pvalues	
	result = pvalue_calculator(N112,N122,N11,N12,theta1 = 0.25,theta2 = 0.5)	
	return(list("nfamily" = nfamily,"nsnp"=nsnp,"pvalue_pr" = result$pvalue_pr,"pvalue_lr" = result$pvalue_lr,"pvalue_lr2"= result$pvalue_lr2,"N11" = N11,"N12" = N12,"N112" = N112,"N122" = N122,"nmissing"= nmissing,"nMedErr" = nMedErr))
}





