waverr <-
function(RawData,Nrepeats){


##PART 1: Create regression formulae##

# Create set of regression formulae (all possible combinations)

Nvariables<-(ncol(RawData))-1

MatrixIntercepts<-matrix(0,nrow=Nvariables,ncol=Nvariables)
MatrixCoefs<-matrix(0,nrow=Nvariables,ncol=Nvariables)
MatrixResiduals<-matrix(NA,nrow=nrow(RawData),ncol=Nvariables^2)
MatrixResidualsStDev<-matrix(0,nrow=Nvariables,ncol=Nvariables)
for (i in 1:Nvariables){
	for (j in 1:Nvariables){
		RegrForm<-lm(RawData[,i+1]~RawData[,j+1])
		MatrixIntercepts[i,j]<-coef(RegrForm)[1]
		MatrixCoefs[i,j]<-coef(RegrForm)[2]
		for (k in 1:length(RegrForm$residuals)){
			MatrixResiduals[k,i+(j-1)*Nvariables]<-RegrForm$residuals[k]
		}	
		MatrixResidualsStDev[i,j]<-sd(RegrForm$residuals)
	}
}

# Calculate log likelihoods

LnLMatrix<-matrix(0,nrow=Nvariables,ncol=Nvariables)
WeightMatrix<-matrix(0,nrow=Nvariables,ncol=Nvariables)


for (i in 1:Nvariables){
	for (j in 1:Nvariables){
		SampleSize<-sum(!is.na(MatrixResiduals[,i+(j-1)*Nvariables]))
		LnL<-(-SampleSize/2)*(log(sum((MatrixResiduals[1:SampleSize,i+(j-1)*Nvariables]^2))/SampleSize))
		LnLMatrix[i,j]<-LnL
	}
	minLnL<-min(LnLMatrix[i,-i])
	maxLnL<-max(LnLMatrix[i,-i])
	for (k in 1:Nvariables){
		Likelihood<-exp(LnLMatrix[i,k])
		SumLikelihoods<-sum(exp(LnLMatrix[i,-i]))
		WeightedLikelihood<-Likelihood/SumLikelihoods		
		WeightMatrix[i,k]<-WeightedLikelihood
	}
}
#################################################################
##PART 2: Reconstruct data using above formulae##

# Figure out which estimations are required for each specimen
print("START reconstruction")
ReconstructedData<-RawData
ReconstructedDataStDev<-matrix(0,nrow=nrow(RawData),ncol=ncol(RawData))
for (i in 1:nrow(RawData)){
	#Work out which variable is present of absent
	ListOfMissingVars<-vector()
	ListOfPresentVars<-vector()
	for (j in 2:ncol(RawData)){
		if (is.na(RawData[i,j])==TRUE) {
			ListOfMissingVars<-c(ListOfMissingVars,j-1)
			}
		else{
			ListOfPresentVars<-c(ListOfPresentVars,j-1)
			}
		}

	if (length(ListOfMissingVars)>0){	
		#Work through each missing variable
		#only work on specimens with missing data
		for (k in 1:length(ListOfMissingVars)){
			VarToPredict<-ListOfMissingVars[k]
			TotalOfEstimations<-0
			TotalOfEstimationsStDev<-0
			TotalWeight<-0
			for (l in 1:length(ListOfPresentVars)){
				PredictorVar<-(ListOfPresentVars[l])
				Estimation<-((MatrixCoefs[VarToPredict,PredictorVar])*RawData[i,PredictorVar+1])+MatrixIntercepts[VarToPredict,PredictorVar]
				EstimationStDev<-MatrixResidualsStDev[VarToPredict,PredictorVar]
				WeightedEstimation<-((WeightMatrix[ListOfMissingVars[k],ListOfPresentVars[l]])+1e-100)*Estimation
				WeightedEstimationStDev<-(WeightMatrix[ListOfMissingVars[k],ListOfPresentVars[l]])*EstimationStDev
				TotalOfEstimations<-TotalOfEstimations+WeightedEstimation
				TotalOfEstimationsStDev<-TotalOfEstimationsStDev+WeightedEstimationStDev
				TotalWeight<-TotalWeight+(WeightMatrix[ListOfMissingVars[k],ListOfPresentVars[l]]+1e-100)
				#print(WeightMatrix[ListOfMissingVars[k],ListOfPresentVars[l]])
			}
			MeanEstimation<-TotalOfEstimations/TotalWeight
			MeanEstimationStDev<-TotalOfEstimationsStDev/TotalWeight
			ReconstructedData[i,ListOfMissingVars[k]+1]<-MeanEstimation
			ReconstructedDataStDev[i,ListOfMissingVars[k]+1]<-MeanEstimationStDev
			#print(paste("Spec.",i,": Var.",VarToPredict,": ",MeanEstimation))
		}
	}
}


# Save files
write.table(ReconstructedData,"ReconstructedData.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(ReconstructedDataStDev,"ReconstructedDataStDev.txt",sep="\t",quote=FALSE,row.names=FALSE)
print("END reconstruction")

print(Sys.time())
##########################################################
##PART 3: Create n datasets resampled from distribution of residuals weighted means##
# Resample the data using paramaters above

print("START resampling reconstructed data")
ResampledData<-array(0, dim=c(nrow(RawData),(ncol(RawData)-1),Nrepeats))
for (i in 1:Nrepeats){
	for (j in 1:nrow(RawData)){
		for (k in 1:(ncol(RawData)-1)){
			ResampledData[j,k,i]<-rnorm(1,ReconstructedData[j,k+1],ReconstructedDataStDev[j,k])
		}
	}
}

VariableStDev<-array(0, dim=c(nrow(RawData),(ncol(RawData)-1)))
for (l in 1:nrow(RawData)){
	for (m in 1:ncol(RawData)-1){
		VariableStDev[l,m]<-sd(ResampledData[l,m,])
	}
}

write.table(ResampledData,"ResampledReconstructions.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(VariableStDev,"ResampledReconstructionsStandardVariation.txt",sep="\t",quote=FALSE,row.names=FALSE)
print("END resampling reconstructed data")

print(Sys.time())
##########################################################






}
