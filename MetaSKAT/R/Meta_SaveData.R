
#
#	Check the Z, and do imputation
#
#
Meta_SKAT_MAIN_Check_Z<-function(Z, n, id_include, SetID, is_dosage=FALSE, impute.method="fixed"){


	m = ncol(Z)
	# OUT

	#############################################
	# Check parameters

	if (class(Z)!= "matrix") stop("Z is not a matrix")
	if (nrow(Z)!=n) stop("Dimensions of y and Z do not match")
 	if(is_dosage ==TRUE){
		impute.method="fixed"
	}

	##############################################
	# Check Missing 

	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 

	###################################################
	# Check missing rates and exclude any SNPs with missing rate > missing_cutoff
	# Also exclude non-polymorphic SNPs
	
	missing.ratio<-rep(0,m)
	for(i in 1:m){
		missing.ratio[i]<-length(which(is.na(Z[,i])))/n
		
		# if missing rate is too high, treat it as unobserved.
		if(missing.ratio[i] > 0.9){
			Z[,i]<-0
		}
	}

	##################################################################
	# doing imputation

	MAF<-colMeans(Z, na.rm = TRUE)/2
	IDX.Err<-which(MAF > 0.5)	
	if(length(IDX.Err) > 0){

		msg<-sprintf("Genotypes of some variants are not the number of minor alleles!")
		#warning(msg,call.=FALSE)
	}
	Z<-SKAT:::Impute(Z,impute.method)

	###########################################
	# Check missing of y and X
	Z<-cbind(Z)
	#id_include1<<-id_include	
	Z.test<-as.matrix(Z[id_include,])

	return(list(Z.test=Z.test, missing_rate=missing.ratio, MAF=MAF , return=0) )

}

Meta_SKAT_SaveData = function(Z, obj.res, SetID, impute.method = "fixed"){

	
	n<-dim(Z)[1]
	m<-dim(Z)[2]
	re<-1

	out.z<-Meta_SKAT_MAIN_Check_Z(Z, n, obj.res$id_include, SetID, impute.method=impute.method)
	
	if(class(obj.res)== "SKAT_NULL_Model_EMMAX"){
		out = Meta_SKAT_SaveData_Kinship(obj.res$res, out.z$Z.test, obj.res$P)
	} else if(obj.res$out_type == "C"){
		out = Meta_SKAT_SaveData_Linear(obj.res$res,out.z$Z.test
			,obj.res$X1, obj.res$s2, obj.res$res.out)
		  
	} else if (obj.res$out_type == "D"){

		out = Meta_SKAT_SaveData_Logistic(obj.res$res, out.z$Z.test
			,obj.res$X1, obj.res$pi_1, obj.res$res.out)
		
	}

	re=list(Score=out$Score, SMat.Summary = out$SMat.Summary, MAF=out.z$MAF, missing_rate=out.z$missing_rate,
	Score.Resampling=out$Score.Resampling )
	
	return(re)

}


Meta_SKAT_SaveData_Linear = function(res, Z, X1, s2, res.out=NULL){

  # get Q
  Q.Temp = t(res)%*%Z

  Q.Temp.Resampling<-NULL
  if(!is.null(res.out)){
 	Q.Temp.Resampling<-t(Z) %*% res.out /s2
  }
  W.1 = t(Z) %*% Z - (t(Z) %*%X1)%*%solve(t(X1)%*%X1)%*% (t(X1) %*% Z ) # t(Z) P0 Z
  MAF = colMeans(Z)/2

  
  #re<-list( Score=Q.Temp[1,] /sqrt(s2), SMat.Summary = W.1, MAF=MAF, Score.Resampling=Q.Temp.Resampling)
  re<-list( Score=Q.Temp[1,] /s2, SMat.Summary = W.1/s2, MAF=MAF, Score.Resampling=Q.Temp.Resampling)  
  return(re)

}

Meta_SKAT_SaveData_Logistic = function(res, Z, X1, pi_1, res.out=NULL){

 
  # Get temp
  Q.Temp = t(res)%*%Z
  Q.Temp.Resampling<-NULL
  if(!is.null(res.out)){
 	Q.Temp.Resampling<-t(Z) %*% res.out
  }

  W.1 = t(Z) %*% (Z * pi_1) - (t(Z * pi_1) %*%X1)%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z * pi_1)) # t(Z) P0 Z
  MAF = colMeans(Z)/2

  re<-list(Score=Q.Temp[1,], SMat.Summary = W.1,  MAF=MAF, Score.Resampling=Q.Temp.Resampling)  
  return(re)
}

Meta_SKAT_SaveData_Kinship = function(res, Z, P1){


  # get Q
  Q.Temp = t(res)%*%Z

  Q.Temp.Resampling<-NULL
  W.1 = t(Z) %*% (P1 %*% Z) # t(Z) P0 Z
  MAF = colMeans(Z)/2

  re<-list( Score=Q.Temp[1,] , SMat.Summary = W.1, MAF=MAF, Score.Resampling=Q.Temp.Resampling)  
  return(re)


}



