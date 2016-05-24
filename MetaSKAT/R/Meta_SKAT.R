
Beta.Weights<-function(MAF,weights.beta, Cutoff=1, Is.MAF=TRUE){

	n<-length(MAF)
	weights<-rep(0,n)
	Sign<-rep(1,n)
	#print(MAF)
	
	IDX1<-which(MAF > 0.5)
	if(length(IDX1) > 0){
		Sign[IDX1]<--1
		MAF[IDX1]<-1-MAF[IDX1]
	}
	 

	
	IDX_0<-union(which(MAF == 0), which(MAF > Cutoff))
	if(length(IDX_0) == n){
		#stop("No polymorphic SNPs")
		weights<-rep(0,n)
	} else if( length(IDX_0) == 0){
		weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
	} else {
		weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
	}

	weights = weights * Sign	
	#if(!Is.MAF){
	#	weights1<<-weights
	#	MAF1<<-MAF
	#} else {
	#	weights2<<-weights
	#	MAF2<<-MAF
	#}
	

	#print(length(IDX_0))
	#print(weights[-IDX_0])
	return(weights)
	
}



Meta_SKAT.Work<-function(re, n.g, combined.weight=TRUE, n1=NULL, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate=FALSE, Group_Idx=NULL, MAF.cutoff=1, Is.MAF=TRUE, missing_cutoff=0.15){


	# optimal =optimal.mod
	if(method=="optimal"){
		method="optimal.mod"
	}


	# Combined MAF
	p<-length(re[[1]]$MAF)
	
	MAF.Combine=0
	MAF.Groups<-list()
	Map.Groups<-rep(0,n.g)


	#re<<-re
	MAF.list<-list()
	for(i in 1:n.g){
	
		MAF.list[[i]]<-re[[i]]$MAF
		
	}
	#MAF.list1<<-MAF.list
	for(i in 1:n.g){
	
		MAF.Combine = MAF.Combine + MAF.list[[i]] * n1[i] / sum(n1)
	}
	
	# If MAF.Combined==0 for all SNP, return p-value 1
	if(sum(MAF.Combine) == 0){
		warning("No polymorphic SNPs!",call.=FALSE)
		return(list(p.value=1, p.value.resampling= NULL, pval.zero.msg=NULL))
	}
	

	# Get MAF.Groups when Group_Idx != NULL
	ID.Groups = unique(Group_Idx)	
	for(j in 1:length(ID.Groups)){
		MAF.Groups[[j]] = 0;
		temp<-which(Group_Idx == ID.Groups[j])
		Map.Groups[temp]<-j
		for(i in temp){
			MAF.Groups[[j]] = MAF.Groups[[j]] + MAF.list[[i]] * n1[i] / sum(n1[temp])
		}
	}


	for(i in 1:n.g){
		if(combined.weight == TRUE){
			weight1<-Beta.Weights(MAF.Combine,weights.beta, MAF.cutoff, Is.MAF=Is.MAF)
		} else {
			j<-Map.Groups[i]
			weight1<-Beta.Weights(MAF.Groups[[j]],weights.beta, MAF.cutoff, Is.MAF=Is.MAF)
		} 
		
		# if missing < missing_cutoff
		
		idx_missing_exclude<-which(re[[i]]$MissingRate > missing_cutoff)
		if(length(idx_missing_exclude) > 0){
			weight1[idx_missing_exclude]<-0
		}
		
		re[[i]]$Score =  re[[i]]$Score * weight1
		re[[i]]$SMat.Summary =  t(t(re[[i]]$SMat.Summary * weight1) * weight1)

		if(!is.null(re[[i]]$Score.Resampling)){
			re[[i]]$Score.Resampling =  re[[i]]$Score.Resampling * weight1
		}
	}

	re.method<-SKAT:::SKAT_Check_Method(method,r.corr)

	if(!is.separate){
		re.score<-Meta_SKAT.Work.OneUnit(re, n.g)
	} else {
		re.score<-Meta_SKAT.Work.Groups(re, n.g, ID.Groups, Group_Idx)
	}


	re<-Met_SKAT_Get_Pvalue(re.score$Score, re.score$SMat.Summary, re.method$r.corr, re.method$method, re.score$Score.Resampling)
	return(re)
	
}


####################################################
#
# Assume SMat and Setinfo are already aligned

MetaSKAT_withlist<-function(SMat.list, Info.list, n.cohort, n.each, combined.weight=TRUE, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate = FALSE, Group_Idx=NULL, MAF.cutoff=1, missing_cutoff=0.15){

	#Info.list1<<-Info.list
	re<-list()
	p<-length(Info.list[[1]]$MAF)
	for(i in 1:n.cohort){
		
		idx_miss<-which(is.na(Info.list[[i]]$MAF))
		if(length(idx_miss) > 0){
			Info.list[[i]]$Score[idx_miss] = 0
			Info.list[[i]]$MAF[idx_miss] = 0
			Info.list[[i]]$MissingRate[idx_miss] = 1
		}

		re1<-list( Score=Info.list[[i]]$Score, SMat.Summary = SMat.list[[i]], MAF=Info.list[[i]]$MAF, MissingRate=Info.list[[i]]$MissingRate)  
		re[[i]]<-re1	
		
	}

	if(is.null(Group_Idx)){
		Group_Idx<-1:n.cohort
	}

	# Use AlleleFreq2 instead of MAF, hence Is.MAF=FALSE
	re = Meta_SKAT.Work(re, n.cohort, combined.weight, n1=n.each, weights.beta=weights.beta, method=method, r.corr=r.corr, is.separate=is.separate, Group_Idx=Group_Idx, 
	MAF.cutoff=MAF.cutoff, Is.MAF=TRUE, missing_cutoff=missing_cutoff )

	return(re)

}

MetaSKAT_MSSD_OneSet<-function(Cohort.Info, SetID, combined.weight=TRUE, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate = FALSE, Group_Idx=NULL, MAF.cutoff=1, missing_cutoff=0.15){

	n.cohort = Cohort.Info$n.cohort
	n.each=rep(0,n.cohort)
	for(i in 1:n.cohort){
		n.each[i]=Cohort.Info$EachInfo[[i]]$header$N.Sample
	}
	
	
	temp<-Get_META_Data_OneSet(Cohort.Info, SetID)
	temp1<-Get_META_Data_OneSet_Align(temp$SMat.list, temp$Info.list, temp$IsExistSNV, n.cohort)

	
	re<-MetaSKAT_withlist(temp1$SMat.list, temp1$Info.list, n.cohort, n.each, combined.weight=combined.weight, 
	weights.beta=weights.beta, method=method, r.corr= r.corr, is.separate = is.separate, Group_Idx=Group_Idx, 
	MAF.cutoff=MAF.cutoff, missing_cutoff=missing_cutoff)
	
	return(re)

}


MetaSKAT_MSSD_ALL<-function(Cohort.Info, ...){

	n.cohort = Cohort.Info$n.cohort
	n.set<-length(Cohort.Info$Set_unique)
	
	pval<-rep(NA,n.set)
	for(i in 1:n.set){
		SetID=Cohort.Info$Set_unique[i]
		out<-try(MetaSKAT_MSSD_OneSet(Cohort.Info,SetID, ... ), silent = TRUE)
		if(class(out)!= "try-error"){
			pval[i]<-MetaSKAT_MSSD_OneSet(Cohort.Info,SetID, ... )$p.value
		}
	}
	
	re<-data.frame(SetID=Cohort.Info$Set_unique, p.value=pval)
	return(re)

}


##################################################################
#
#	Genotype matrix Z should be matched with y and X
#

# change from here
MetaSKAT_wZ<-function(Z, obj, combined.weight=TRUE, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate = FALSE, Group_Idx=NULL, impute.method="fixed", impute.estimate.maf=1, missing_cutoff=0.15){


	if(is.matrix(Z)!= TRUE){
		stop("ERROR: Z is not a matrix!")
	}
	if(class(obj)!= "META_NULL_Model" && class(obj)!= "META_NULL_Model_EmmaX"){
		stop("ERROR: obj class is not either META_NULL_Model or META_NULL_Model_EmmaX!")
	}
	
	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 
	
	nSNP<-ncol(Z)
	n<-nrow(Z)
	n.g<-obj$n.g
	
	# Calculate missing rate : Since Z is imputed before calling Meta_SKAT_SaveData when impute.estimate.maf==1, 
	# missing rate should be calculated before calling this function.
	
	missing.ratio.Z = matrix(rep(0, nSNP*n.g), nrow=n.g)
	for(i in 1:n.g){
		ID<-obj$ID[[i]]
		Z1<-as.matrix(Z[ID,])
		n1<-nrow(Z1)
		for(j in 1:nSNP){
			missing.ratio.Z[i,j]<-length(which(is.na(Z1[,j])))/n1
		}
	}
	
	Is.impute.cohortwise= FALSE
	if(length(IDX_MISS) > 0){

		msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )

		warning(msg,call.=FALSE)
		
		if(impute.estimate.maf==1){
			Is.impute.cohortwise = TRUE
		
		} else if(impute.estimate.maf==2){
			Z<-SKAT:::Impute(Z,impute.method=impute.method)
		} else {
			stop("ERROR: impute.estimate.mat is wrong! it should be either 1 or 2")
		}
	} 
	
	
	re1<-list()
	for(i in 1:n.g){

		ID<-obj$ID[[i]]
		Z1<-as.matrix(Z[ID,])
		re1[[i]]<-Meta_SKAT_SaveData(Z1, obj$out[[i]], SetID=NULL, impute.method = impute.method)
		
		re1[[i]]$MissingRate = missing.ratio.Z[i,]
	}

	if(is.null(Group_Idx)){
		Group_Idx<-1:n.g
	}

	re = Meta_SKAT.Work(re1, n.g, combined.weight, n1=obj$n.each, weights.beta=weights.beta, method=method, r.corr=r.corr,is.separate=is.separate
	, Group_Idx=Group_Idx, missing_cutoff=missing_cutoff)

	return(re)

}





MetaSKAT_wZ_OLD<-function(Z, obj, combined.weight=TRUE, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate = FALSE, Group_Idx=NULL, impute.method="fixed", impute.estimate.maf=1, missing_cutoff=0.15){


	if(is.matrix(Z)!= TRUE){
		stop("ERROR: Z is not a matrix!")
	}

	if(class(obj)!= "META_NULL_Model" && class(obj)!= "META_NULL_Model_EmmaX"){
		stop("ERROR: obj class is not either META_NULL_Model or META_NULL_Model_EmmaX!")
	}
	
	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 
	
	nSNP<-ncol(Z)
	n<-nrow(Z)
	Is.impute.cohortwise= FALSE
	missing.ratio.Z = rep(0, nSNP)
	for(i in 1:nSNP){
		missing.ratio.Z[i]<-length(which(is.na(Z[,i])))/n
	}
	
	if(length(IDX_MISS) > 0){

		msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )

		warning(msg,call.=FALSE)
		
		if(impute.estimate.maf==1){
			Is.impute.cohortwise = TRUE
		
		} else if(impute.estimate.maf==1){
			Z<-SKAT:::Impute(Z,impute.method=impute.method)
		} else {
			stop("ERROR: impute.estimate.mat is wrong! it should be either 1 or 2")
		}
	} 
	
	n.g<-obj$n.g
	re1<-list()
	for(i in 1:n.g){

		ID<-obj$ID[[i]]
		if(Is.impute.cohortwise ){
			Z1<-as.matrix(Z[ID,])
		} else {
			Z1<-SKAT:::Impute(as.matrix(Z[ID,]),impute.method=impute.method)
		}
		res<-obj$out[[i]]$res
		res.out<-obj$out[[i]]$res.out
		X1<-obj$out[[i]]$X1
		

		if(obj$out_type=="C"){
			s2<-obj$out[[i]]$s2
			re1[[i]]<-Meta_SKAT_SaveData_Linear(res,Z1 , X1, s2, res.out)
		} else if (obj$out_type=="D"){
			pi_1<-obj$out[[i]]$pi_1
			re1[[i]]<-Meta_SKAT_SaveData_Logistic(res,Z1 , X1, pi_1, res.out)
		} else if (obj$out_type=="K"){
			re1[[i]]<-Meta_SKAT_SaveData_Kinship(res,Z1 , obj$out[[i]]$P)
		} else {
			stop("ERROR: out_type is wrong!")
		}
		
	}

	if(is.null(Group_Idx)){
		Group_Idx<-1:n.g
	}

	re = Meta_SKAT.Work(re1, n.g, combined.weight, n1=obj$n.each, weights.beta=weights.beta, method=method, r.corr=r.corr,is.separate=is.separate
	, Group_Idx=Group_Idx, missing_cutoff=missing_cutoff)

	return(re)

}



