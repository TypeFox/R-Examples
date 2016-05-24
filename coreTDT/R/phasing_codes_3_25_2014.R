##NOTE: This is for simulations. This program cannot deal with sequencing errors, thus cannot apply to real data analysis


##parent geneotype [,,1] for father [,,2] for mother

##deal with parent-child pair
PairPhase = function(paternal.genotype,child.genotype){
	paternal.return = 0
	child.return = 0
	if(all(paternal.genotype ==0)){
		return(c(paternal.return, child.return))
	}else{
		p.c.prod = paternal.genotype * child.genotype
		if(any(p.c.prod !=0)){
			child.return = 1
		}
		p.loci = which(paternal.genotype !=0)
		paternal.return = (0 %in% p.c.prod[p.loci]) + (any(p.c.prod[p.loci] !=0))
		
		return(c(paternal.return, child.return))
	}
	
}


##deal with a family, considering homozygous mutations
FamilyPhase = function(parent.genotype,child.genotype){
	nsnv = length(child.genotype)
	##check if 2 exist in parents
	##both parents are homozygous
	if((2 %in% parent.genotype[,2]) & (2 %in% parent.genotype[,1])){
		maternal.return = 2
		paternal.return = 2	
		child.return = 2
	}else if(2 %in% parent.genotype[,1]){
		##only father has homozygous mutatios
		paternal.return = 2
		homo.loci = which(parent.genotype[,1] == 2)
		tmp.vec = rep(0,nsnv)
		tmp.vec[homo.loci] = 1
		child.genotype.tmp = child.genotype - tmp.vec
		##check if mistakes occurs
		if(all(child.genotype.tmp[homo.loci] %in% c(0,1))==F){
			stop("mistakes in substracting homozygous mutatios")
		}
		
		tmp.result = PairPhase(parent.genotype[,2],child.genotype.tmp)
		maternal.return = tmp.result[1]
		child.return = tmp.result[2]+1
	}else if(2 %in% parent.genotype[,2]){
		##only mother has homozygous mutatios
		maternal.return = 2
		homo.loci = which(parent.genotype[,2] == 2)
		tmp.vec = rep(0,nsnv)
		tmp.vec[homo.loci] = 1
		child.genotype.tmp = child.genotype - tmp.vec
		##check if mistakes occurs
		if(all(child.genotype.tmp[homo.loci] %in% c(0,1))==F){
			stop("mistakes in substracting homozygous mutatios")
		}
		
		tmp.result = PairPhase(parent.genotype[,1],child.genotype.tmp)
		paternal.return = tmp.result[1]
		child.return = tmp.result[2]+1
	}else{
		#no homozygous mutations
		tmp.result1 = PairPhase(parent.genotype[,1],child.genotype)
		tmp.result2 = PairPhase(parent.genotype[,2],child.genotype)
		paternal.return = tmp.result1[1]
		maternal.return = tmp.result2[1]
		child.return = tmp.result1[2] + tmp.result2[2]
	}
	
	return(c(paternal.return, maternal.return, child.return))
	
}

## including the situation with shared mutations

##if any shared loci are heterozygous in all family memebers, then they will not be able to produce any information
FamilyPhaseII = function(parent.genotype,child.genotype){
	
	
	rm.loci =  which(parent.genotype[,1] ==1 &parent.genotype[,2] ==1 )
		
	if(length(rm.loci) !=0){
		if( (length(rm.loci) ==1) & (sum(c(parent.genotype[-rm.loci,])) == 0) ){		
			result.return = c(parent.genotype[rm.loci,], child.genotype[rm.loci])
		}else{
			result.return = FamilyPhase(parent.genotype[-rm.loci,], child.genotype[-rm.loci])
		}	
	}else{
		result.return = FamilyPhase(parent.genotype, child.genotype)
	}
		
	return(result.return)	
}


##compare to FamilyPhaseII, this one directly set families with shared heterzygous mutations as missing data
FamilyPhaseIII = function(parent.genotype,child.genotype){
	
	
	rm.loci =  which(parent.genotype[,1] ==1 &parent.genotype[,2] ==1 )
		
	##only one mutation in the family and it is heterozygous shared by parents	
	if(length(rm.loci) !=0){	
		if( (length(rm.loci) ==1) & (sum(c(parent.genotype[-rm.loci,])) == 0) ){		
			result.return = c(parent.genotype[rm.loci,], child.genotype[rm.loci])
		}else{
			##other case just set the family as missing data
			result.return = rep(0,3)
		}	
	}else{
		result.return = FamilyPhase(parent.genotype, child.genotype)
	}
		
	return(result.return)	
}












