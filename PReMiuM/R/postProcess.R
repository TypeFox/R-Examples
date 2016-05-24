# (C) Copyright David Hastie, Silvia Liverani and Aurore J. Lavigne, 2012-2014.

# PReMiuM++ is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.

# PReMiuM++ is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with PReMiuM++ in the documentation directory. If not, see
# <http://www.gnu.org/licenses/>.

# The external linear algebra library Eigen, parts of which are included  in the
# lib directory is released under the LGPL3+ licence. See comments in file headers
# for details.

# The Boost C++ header library, parts of which are included in the  lib directory
# is released under the Boost Software Licence, Version 1.0, a copy  of which is
# included in the documentation directory.

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

profRegr<-function(covNames, fixedEffectsNames, outcome="outcome", outcomeT=NA, data, output="output", hyper, predict, predictType="RaoBlackwell", nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1, nClusInit, seed, yModel="Bernoulli", xModel="Discrete", sampler="SliceDependent", alpha=-2, dPitmanYor=0, excludeY=FALSE, extraYVar=FALSE, varSelectType="None", entropy,reportBurnIn=FALSE, run=TRUE, discreteCovs, continuousCovs, whichLabelSwitch="123", includeCAR=FALSE, neighboursFile="Neighbours.txt",weibullFixedShape=TRUE, useNormInvWishPrior=FALSE){

	# suppress scientific notation
	options(scipen=999)

	if (xModel=="Mixed"){
		covNames <- c(discreteCovs, continuousCovs)
		nDiscreteCovs <- length(discreteCovs)
		nContinuousCovs <- length(continuousCovs)
	}
	nCovariates<-length(covNames)
	
	if (!is.data.frame(data)) stop("Input data must be a data.frame with outcome, covariates and fixed effect names as column names.")

	if (extraYVar==TRUE&(yModel=="Categorical"||yModel=="Normal"||yModel=="Survival")) stop("Option extraYVar is only available for Bernoulli, Binomial and Poisson response.")

	if (includeCAR==TRUE&(yModel=="Categorical"||yModel=="Survival"||yModel=="Bernoulli"||yModel=="Binomial")) stop("Option includeCAR is only available for Poisson and Normal response.")

	if (predictType=="random"){
		if (yModel!="Normal") stop("The option of random predictions is only available for yModel=Normal.")
	}

	if (useNormInvWishPrior==TRUE && !varSelectType=="None") stop("Variable selection is not available for Normal-inverse-Wishart prior for Normal covariates.")

	# open file to write output
	fileName<-paste(output,"_input.txt",sep="")
	# make big data matrix with outcome, covariates and fixed effects	
	# outcome
	# create outcome if excludeY=TRUE and outcome not provided
	if (length(which(colnames(data)==outcome))<1&&excludeY==TRUE) {
		dataMatrix<-rep(0,dim(data)[1])
		yModel="Bernoulli" 
	} else {
		dataMatrix<-data[,which(colnames(data)==outcome)]
	}

	if (sum(is.na(dataMatrix))>0) stop("ERROR: the outcome cannot have missing values. Use the profiles with missing outcome for predictions.")

	if (length(dataMatrix)==0) stop("ERROR: the outcome provided is not one of the names of the variables in the dataframe.")

	# recode outcome covariates
	if (yModel=="Categorical"||yModel=="Bernoulli"){
		outcomeY<-dataMatrix
		outcomeFactor<-as.factor(outcomeY)
		yLevels<-length(levels(outcomeFactor))
		if (yModel=="Bernoulli"&yLevels>2) stop("The number of levels of the outcome is greater than 2, which is not allowed for Bernoulli outcome. You might want to set yModel to be Categorical.") 
		if (yModel=="Categorical") write(yLevels,fileName,append=T,ncolumns=length(yLevels))
		if (is.numeric(outcomeY)){
			if (!(min(outcomeY)==0&&max(outcomeY)==(yLevels-1)&&sum(!is.wholenumber(outcomeY))==0)) {
				print("Recoding of the outcome as follows")
				print(paste("Replacing level ",levels(outcomeFactor)," with ",c(0:(yLevels-1)),sep=""))
				levels(outcomeFactor)<-c(0:(yLevels-1))
				dataMatrix<-outcomeFactor
			}
		} else {
			print("Recoding of the outcome as follows")
			tmpLevels<-levels(outcomeFactor)
			print(paste("Replacing level ",levels(outcomeFactor)," with ",c(0:(yLevels-1)),sep=""))
			levels(outcomeFactor)<-c(0:(yLevels-1))
			dataMatrix<-outcomeFactor
		}
		if (yModel=="Bernoulli") yLevels <- 1
	} else {
		yLevels <- 1
	}

	# covariates
	covIndeces<-vector()
	for (i in 1:nCovariates){
		tmpIndex<-which(colnames(data)==covNames[i])
		if (length(tmpIndex)==0) stop("ERROR: covariate names in data.frame provided do not correspond to list of covariates for profile regression")
		covIndeces<-append(covIndeces,tmpIndex)
	}

	covariates<-data[,covIndeces]
	dataMatrix<-cbind(dataMatrix,covariates)

	# recode covariate levels
	if (xModel=="Discrete"){
		xLevels<-vector()
		for (k in 1:nCovariates){
			tmpCov<-dataMatrix[,(1+k)]
			xLevels[k]<-length(levels(as.factor(tmpCov)))	
			if (!(min(tmpCov,na.rm=TRUE)==0&&max(tmpCov,na.rm=TRUE)==(xLevels[k]-1)&&sum(!is.wholenumber(tmpCov[!is.na(tmpCov)]))==0)) {
				print(paste("Recoding of covariate ",colnames(dataMatrix)[k+1]," as follows",sep=""))
				tmpCovFactor<-as.factor(tmpCov)
				tmpLevels<-levels(tmpCovFactor)
				print(paste("Replacing level ",levels(tmpCovFactor)," with ",c(0:(xLevels[k]-1)),sep=""))
				levels(tmpCovFactor)<-c(0:(xLevels[k]-1))	
				dataMatrix[,(1+k)]<-tmpCovFactor
				dataMatrix[,(1+k)]<-as.numeric(levels(dataMatrix[,(1+k)]))[as.integer(dataMatrix[,(1+k)])]
			}
		}
	} else 	if (xModel=="Mixed"){
		xLevels<-vector()
		for (k in 1:nDiscreteCovs){
			tmpCov<-dataMatrix[,(1+k)]
			xLevels[k]<-length(levels(as.factor(tmpCov)))	
			if (!(min(tmpCov,na.rm=TRUE)==0&&max(tmpCov,na.rm=TRUE)==(xLevels[k]-1)&&sum(!is.wholenumber(tmpCov[!is.na(tmpCov)]))==0)) {
				print(paste("Recoding of covariate number ",colnames(dataMatrix)[k+1]," as follows",sep=""))
				tmpCovFactor<-as.factor(tmpCov)
				tmpLevels<-levels(tmpCovFactor)
				print(paste("Replacing level ",levels(tmpCovFactor)," with ",c(0:(xLevels[k]-1)),sep=""))
				levels(tmpCovFactor)<-c(0:(xLevels[k]-1))	
				dataMatrix[,(1+k)]<-tmpCovFactor
				dataMatrix[,(1+k)]<-as.numeric(levels(dataMatrix[,(1+k)]))[as.integer(dataMatrix[,(1+k)])]
			}
		}
	} else {
		xLevels <- 1
	}

	for (k in 1:nCovariates){
		missingX<-is.na(dataMatrix[,(k+1)])
		nMissingX<-sum(missingX)
		if (nMissingX>0) {
			dataMatrix[missingX,(k+1)]<- rep(-999,nMissingX)
		}
	}

	# fixed effects
	if (!missing(fixedEffectsNames)) {
		nFixedEffects<-length(fixedEffectsNames)
		FEIndeces<-vector()
		for (i in 1:nFixedEffects){
			tmpIndex<-which(colnames(data)==fixedEffectsNames[i])
			if (length(tmpIndex)==0) stop("ERROR: fixed effects names in data.frame provided do not correspond to list of fixed effects for profile regression")
			FEIndeces<-append(FEIndeces,tmpIndex)
		}
		fixedEffects<-data[,FEIndeces]
		if (sum(is.na(fixedEffects))>0) stop("ERROR: fixed effects cannot have missing values. Use an imputation method before using profRegr().")
		dataMatrix<-cbind(dataMatrix,fixedEffects)
		for (i in dim(fixedEffects)[2]){
			if (class(fixedEffects[,i])=="character") stop("ERROR: fixed effects must be of class numeric. See help pages.") 
		}
	} else {
		nFixedEffects<-0
		#if (yModel=="Survival") stop("ERROR: For the current implementation of Survival outcome the fixed effects must be provided. ")
	}

	#  extra outcome data
	if (yModel=="Poisson"||yModel=="Binomial"||yModel=="Survival") {
		if(is.na(outcomeT)){
			stop ("It is required to set outcomeT for Poisson (offset), Binomial (number of trials) or Survival (censoring) outcome.")
		} else {
			indexOutcomeT <- which(colnames(data)==outcomeT)
			if (yModel=="Survival") {
				if (!all(is.element(names(table(data[indexOutcomeT])),c("0","1")))) stop("If yModel is Survival, variable outcomeT corresponds to the censoring of the outcome. Therefore, it must be coded 0 and 1.")

			}
			dataMatrix <- cbind(dataMatrix,data[indexOutcomeT])
		}
	} else {
		if(!is.na(outcomeT)) stop ("It is only required to set outcomeT for Poisson, Binomial and Survival outcome.")
	}		

	# print number of subjects
	nSubjects <- dim(dataMatrix)[1]
	write(as.character(nSubjects), fileName,ncolumns=1)
	# print number of covariates and their names
	write(as.character(nCovariates),fileName,append=T,ncolumns=1)
	if (xModel=="Mixed"){
		write(as.character(nDiscreteCovs),fileName,append=T,ncolumns=1)
		write(as.character(nContinuousCovs),fileName,append=T,ncolumns=1)
	}
	write(t(covNames), fileName,append=T,ncolumns=1)
	# print number of fixed effects and their names
	write(nFixedEffects, fileName,append=T,ncolumns=1)
	if (nFixedEffects>0){
		write(t(fixedEffectsNames), fileName,append=T,ncolumns=1)
	}
	if (yModel=="Categorical") write(yLevels,fileName,append=T,ncolumns=1) 
	if (xModel=="Discrete"||xModel=="Mixed"){
		write(xLevels,fileName,append=T,ncolumns=length(xLevels))
	}

	# write prediction file
#	if (!missing(predict)&(includeCAR)) stop ("Predictions are not available with spatial effect.")
	if (!missing(predict)) {
		nPreds<-dim(predict)[1]
		write(nPreds, paste(output,"_predict.txt",sep=""),ncolumns=1)
		for (k in 1:nPreds){
			missingX<-is.na(predict[k,])
			nMissingX<-sum(missingX)
			if (nMissingX>0) {
				predict[k,missingX]<- rep(-999,nMissingX)
			}
		}		
		write(t(as.matrix(predict[,c(covNames)])), paste(output,"_predict.txt",sep=""),append=T,ncolumns=length(covNames))
		if (length(intersect(outcome,names(predict)))>0) {
			if (nFixedEffects>0 ) {
			if (length(intersect(fixedEffectsNames,names(predict)))==length(fixedEffectsNames)) {
				write(nPreds, paste(output,"_predictFull.txt",sep=""),ncolumns=1)
				if (yModel=="Poisson" || yModel=="Binomial"||yModel=="Survival"){
					write(t(as.matrix(predict[,c(outcome,fixedEffectsNames,outcomeT)])), paste(output,"_predictFull.txt",sep=""),append=T,ncolumns=(2+length(fixedEffectsNames)))
				} else {
					write(t(as.matrix(predict[,c(outcome,fixedEffectsNames)])), paste(output,"_predictFull.txt",sep=""),append=T,ncolumns=(1+length(fixedEffectsNames)))				
				}
				fullPredictFile<-TRUE
			} else {
				write(nPreds, paste(output,"_predictFull.txt",sep=""),ncolumns=1)
				predictFixedEffectsNA<-cbind(as.matrix(predict[,c(outcome)]),matrix(-999,ncol=nFixedEffects,nrow=nPreds))
				write(t(predictFixedEffectsNA), paste(output,"_predictFull.txt",sep=""),append=T,ncolumns=(1+length(fixedEffectsNames)))
				fullPredictFile<-TRUE
			}
			} else {

				write(nPreds, paste(output,"_predictFull.txt",sep=""),ncolumns=1)
				if (yModel=="Poisson" || yModel=="Binomial" ||yModel == "Survival"){
					predictFixedEffectsNA<-as.matrix(predict[,c(outcome,outcomeT)])
					write(t(predictFixedEffectsNA), paste(output,"_predictFull.txt",sep=""),append=T,ncolumns=2)
				} else {
					predictFixedEffectsNA<-as.matrix(predict[,c(outcome)])
					write(t(predictFixedEffectsNA), paste(output,"_predictFull.txt",sep=""),append=T,ncolumns=1)
				}
				fullPredictFile<-TRUE
			}

		} else {
			fullPredictFile<-FALSE
		}
	} else {
		nPreds<-0
		fullPredictFile<-FALSE
	}

	write(t(dataMatrix), fileName,append=T,ncolumns=dim(dataMatrix)[2])

	# other checks to ensure that there are no errors when calling the program
	if (xModel!="Discrete"&xModel!="Normal"&xModel!="Mixed") stop("This xModel is not defined.")
	if (yModel!="Poisson"&yModel!="Binomial"&yModel!="Bernoulli"&yModel!="Normal"&yModel!="Categorical"&yModel!="Survival") stop("This yModel is not defined.")

	# conditions for alpha and dPitmanYor parameters	
	# checks that dPitmanYor is in the correct interval
	if (dPitmanYor<0&dPitmanYor>=1) stop("dPitmanYor must belongto the interval [0,1).")
	# if alpha < 0 then we have a Dirichlet prior, so cannot use dPitmanYor
	# if the user asks for alpha<0 and dPitmanYor>0 we stop and give an error message to ensure the user is aware
	if (alpha<=-1&&dPitmanYor>0) stop("Setting alpha < 0 (the default is -2) and dPitmanYor > 0 is not allowed.")
	# the third label switching move is only for Dirichlet process prior, not for Pitman-Yor process prior 
	if (dPitmanYor>0) whichLabelSwitch<-"12"
	if ((alpha>=-1)&(alpha < -dPitmanYor)) stop("It is a condition of the Pitman-Yor process prior that alpha > - dPitmanYor.")
	if (dPitmanYor>0&sampler!="Truncated") stop("The Pitman-Yor process prior has been implemented to use with the Truncated sampler only.")
	#if (dPitmanYor>0&sampler=="Truncated") print("Note that for the Pitman-Yor process prior there might be en error when using the Truncated sampler. This is due to the way that the bound on the number of clusters is computed.")

	#check entries for spatial CAR term
	if (includeCAR&file.exists(neighboursFile)==FALSE) stop("You must enter a valid file for neighbourhood structure.") 
 
	inputString<-paste("PReMiuM --input=",fileName," --output=",output," --xModel=",xModel," --yModel=",yModel," --varSelect=",varSelectType," --whichLabelSwitch=",whichLabelSwitch," --predType=",predictType,sep="")

	# create hyperparameters file
	if (!missing(hyper)) {
		hyperFile <-paste(output,"_hyper.txt",sep="")
		if (file.exists(hyperFile)) file.remove(hyperFile)
		if (!is.null(hyper$shapeAlpha)){
			write(paste("shapeAlpha=",hyper$shapeAlpha,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$rateAlpha)){
			write(paste("rateAlpha=",hyper$rateAlpha,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$aPhi)){
			write(paste("aPhi=",paste(hyper$aPhi,collapse=" ")," ",sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$mu0)){
			write(paste("mu0=",paste(hyper$mu0,collapse=" ")," ",sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$Tau0)){
			write(paste("Tau0=",paste(t(hyper$Tau0),collapse=" ")," ",sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$R0)){
			write(paste("R0=",paste(t(hyper$R0),collapse=" ")," ",sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$kappa0)){
			write(paste("kappa0=",hyper$kappa0,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$nu0)){
		  write(paste("nu0=",hyper$nu0,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$muTheta)){
			write(paste("muTheta=",hyper$muTheta,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$sigmaTheta)){
			write(paste("sigmaTheta=",hyper$sigmaTheta,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$dofTheta)){
			write(paste("dofTheta=",hyper$dofTheta,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$muBeta)){
			write(paste("muBeta=",hyper$muBeta,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$sigmaBeta)){
			write(paste("sigmaBeta=",hyper$sigmaBeta,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$dofBeta)){
			write(paste("dofBeta=",hyper$dofBeta,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$shapeTauEpsilon)){
			write(paste("shapeTauEpsilon=",hyper$shapeTauEpsilon,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$rateTauEpsilon)){
			write(paste("rateTauEpsilon=",hyper$rateTauEpsilon,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$aRho)){
			write(paste("aRho=",hyper$aRho,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$bRho)){
			write(paste("bRho=",hyper$bRho,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$atomRho)){
			if (hyper$atomRho<=0 || hyper$atomRho >1) stop("Hyperparameter atomRho must be in (0,1]. See ?setHyperparams for help.")
			write(paste("atomRho=",hyper$atomRho,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$shapeSigmaSqY)){
			write(paste("shapeSigmaSqY=",hyper$shapeSigmaSqY,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$scaleSigmaSqY)){
			write(paste("scaleSigmaSqY=",hyper$scaleSigmaSqY,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$rSlice)){
			write(paste("rSlice=",hyper$rSlice,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$truncationEps)){
			write(paste("truncationEps=",hyper$truncationEps,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$shapeTauCAR)){
		  write(paste("shapeTauCAR=",hyper$shapeTauCAR,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$rateTauCAR)){
		  write(paste("rateTauCAR=",hyper$rateTauCAR,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$shapeNu)){
			write(paste("shapeNu=",hyper$shapeNu,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$scaleNu)){
			write(paste("scaleNu=",hyper$scaleNu,sep=""),hyperFile,append=T)
		}
		if (!is.null(hyper$initAlloc)){
			write(paste("initAlloc=",paste(hyper$initAlloc,collapse=" ")," ",sep=""),hyperFile,append=T)
			if(sum(is.wholenumber(hyper$initAlloc))!=nSubjects) stop("The vector initAlloc has not been initialised properly: the vector must contain integers and have as many elements as there are observations in the dataset (nSubjects).")
			nClusInit<-max(hyper$initAlloc)

				
		}
	}

	if (weibullFixedShape) inputString<-paste(inputString," --weibullFixedShape",sep="")
	if (reportBurnIn) inputString<-paste(inputString," --reportBurnIn",sep="")
	if (!missing(alpha)) inputString<-paste(inputString," --alpha=",alpha,sep="")
	if (!missing(dPitmanYor)) inputString<-paste(inputString," --dPitmanYor=",dPitmanYor,sep="")
	if (!missing(sampler)) inputString<-paste(inputString," --sampler=",sampler,sep="")
	if (!missing(hyper)) inputString<-paste(inputString," --hyper=",hyperFile,sep="")
	if (!missing(predict)) inputString<-paste(inputString," --predict=",paste(output,"_predict.txt",sep=""),sep="")
	if (!missing(nSweeps)) inputString<-paste(inputString," --nSweeps=",nSweeps,sep="")
	if (!missing(nBurn)) inputString<-paste(inputString," --nBurn=",nBurn,sep="")
	if (!missing(nProgress)) inputString<-paste(inputString," --nProgress=",nProgress,sep="")
	if (!missing(nFilter)) inputString<-paste(inputString," --nFilter=",nFilter,sep="")
	if (!missing(nClusInit)) inputString<-paste(inputString," --nClusInit=",nClusInit,sep="")
	if (!missing(seed)) inputString<-paste(inputString," --seed=",seed,sep="")
	if (excludeY) inputString<-paste(inputString," --excludeY",sep="")
	if (extraYVar) inputString<-paste(inputString," --extraYVar",sep="")
	if (!missing(entropy)) inputString<-paste(inputString," --entropy",sep="")
	if (includeCAR) inputString<-paste(inputString," --includeCAR", " --neighbours=", neighboursFile ,sep="")
	if (useNormInvWishPrior) inputString<-paste(inputString," --useNormInvWishPrior", sep="")

	if (run) .Call('profRegr', inputString, PACKAGE = 'PReMiuM')


	# define directory path and fileStem
	outputSplit <- strsplit(output,split="/")
	fileStem <- tail(outputSplit[[1]],1)
	if (length(outputSplit[[1]])>1) {
		directoryPath <- paste(head(outputSplit[[1]],-1),collapse="/")
	} else {
		directoryPath <-"."
	}
	
	# other re-writes for function return
	# var select and var select type
	if (varSelectType=="None") {
		varSelect <- FALSE
		varSelType <- NULL
	} else {
		varSelect <- TRUE
		varSelType <- varSelectType
	}
	# covariate matrix
	xMat <- dataMatrix[,2:(nCovariates+1)]
	# outcome and fixed effect matrix
	yMat <- NULL
	wMat <- NULL
	# the code requires these matrices whether excludeY is TRUE or FALSE
	#if(includeResponse){
		yMat<-matrix(dataMatrix[,1],ncol=1)
		if(yModel=='Poisson'){
			offset<-dataMatrix[,ncol(dataMatrix)]
			yMat<-cbind(yMat,offset)
		}else if(yModel=='Binomial'){
			nTrials<-dataMatrix[,ncol(dataMatrix)]
			yMat<-cbind(yMat,nTrials)
		}else if(yModel=='Survival'){
			censoring<-dataMatrix[,ncol(dataMatrix)]
			yMat<-cbind(yMat,censoring)
		}
		if(nFixedEffects>0){
			wMat<-dataMatrix[,(2+nCovariates):(1+nCovariates+nFixedEffects)]
		}
	#}
	# include response
	if (excludeY) {
		includeResponse <- FALSE
		yModel <- NULL
	} else {
		includeResponse <- TRUE
	}

	if (xModel=="Mixed") {
		discreteTmp<-discreteCovs
		contTmp<-continuousCovs
	} else {
		discreteTmp<-NA
		contTmp<-NA
	}

	if (xModel=="Discrete"){
		useNIWP<-NA
	}else{
		useNIWP<-useNormInvWishPrior
	}


	return(list("directoryPath"=directoryPath,
		"fileStem"=fileStem,
		"inputFileName"=fileName,
		"nSweeps"=nSweeps,
		"nBurn"=nBurn,
		"reportBurnIn"=reportBurnIn,
		"nFilter"=nFilter,
		"nProgress"=nProgress,
		"nSubjects"=nSubjects,
		"nPredictSubjects"=nPreds,
		"fullPredictFile"=fullPredictFile,
		"predictType"=predictType,
		"alpha"=alpha,
		"dPitmanYor"=dPitmanYor,
		"covNames"=covNames,
		"discreteCovs"=discreteTmp,
		"continuousCovs"=contTmp,
		"xModel"=xModel,
		"includeResponse"=includeResponse,
		"whichLabelSwitch"=whichLabelSwitch,
		"yModel"=yModel,
		"varSelect"=varSelect,
		"varSelectType"=varSelType,
		"nCovariates"=nCovariates,
		"nDiscreteCovs"=ifelse(xModel=="Mixed",nDiscreteCovs,NA),
		"nContinuousCovs"=ifelse(xModel=="Mixed",nContinuousCovs,NA),
		"nFixedEffects"=nFixedEffects,
		"nCategoriesY"=yLevels,
		"nCategories"=xLevels,
		"extraYVar"=extraYVar,
		"includeCAR"=includeCAR,
		"weibullFixedShape"=weibullFixedShape,
		"useNormInvWishPrior"=useNIWP,
		"xMat"=xMat,"yMat"=yMat,"wMat"=wMat))
}



# Function to take the output from the C++ run and return an average dissimilarity
# matrix
calcDissimilarityMatrix<-function(runInfoObj,onlyLS=FALSE){

	directoryPath=NULL
	fileStem=NULL
	reportBurnIn=NULL
	nSweeps=NULL
	nFilter=NULL
	nSubjects=NULL
	nPredictSubjects=NULL
	nBurn=NULL

   for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])

   fileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
 	   
   if (reportBurnIn) {
	recordedNBurn<-nBurn
    } else {
	recordedNBurn<-1
    }

   # Call the C++ to compute the dissimilarity matrix
   disSimList<-.Call('calcDisSimMat',fileName,nSweeps,recordedNBurn,nFilter,nSubjects,
                       nPredictSubjects, onlyLS, PACKAGE = 'PReMiuM')

	if (onlyLS){
		lsOptSweep<-disSimList$lsOptSweep
		disSimMatPred<-NULL              
		disSimObj<-list('disSimRunInfoObj'=runInfoObj,'disSimMat'=NA,
			'disSimMatPred'=NA,'lsOptSweep'=lsOptSweep,'onlyLS'=onlyLS)              
	} else {
		disSimMat<-disSimList$disSimMat
		lsOptSweep<-disSimList$lsOptSweep
		disSimMatPred<-NULL              
		if(nPredictSubjects>0){
			disSimMatPred<-disSimMat[(1+(nSubjects*(nSubjects-1)/2)):length(disSimMat)]
			disSimMat<-disSimMat[1:(nSubjects*(nSubjects-1)/2)]
		}   
		disSimObj<-list('disSimRunInfoObj'=runInfoObj,'disSimMat'=disSimMat,
			'disSimMatPred'=disSimMatPred,'lsOptSweep'=lsOptSweep,'onlyLS'=onlyLS)              
	}

	return(disSimObj)
}

# Given a dissimilarity matrix (or list of dissimilarity matrices)
# run partitioning around medoids clustering
calcOptimalClustering<-function(disSimObj,maxNClusters=NULL,useLS=F){
   
	disSimRunInfoObj=NULL
	directoryPath=NULL
	fileStem=NULL
	lsOptSweep=NULL
	nSubjects=NULL
	nPredictSubjects=NULL
	reportBurnIn=NULL
	nBurn=NULL
	nFilter=NULL
	nSweeps=NULL
	onlyLS=NULL

   for (i in 1:length(disSimObj)) assign(names(disSimObj)[i],disSimObj[[i]])
   for (i in 1:length(disSimRunInfoObj)) assign(names(disSimRunInfoObj)[i],disSimRunInfoObj[[i]])

	if (onlyLS==TRUE) useLS <- TRUE

   if(useLS){
      # maniupulation for least squares method, but computation has been done in previous function
      zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
      zFile<-file(zFileName,open="r")
            
      optZ<-scan(zFile,what=integer(),skip=lsOptSweep-1,n=nSubjects+nPredictSubjects,quiet=T)
      optZFit<-optZ[1:nSubjects]
      if(nPredictSubjects>0){
         optZPredict<-optZ[(nSubjects+1):(nSubjects+nPredictSubjects)]
      }
      uZFit<-unique(optZFit)
      chosenNClusters<-length(unique(uZFit))
      clustVec<-match(optZFit,uZFit)
      clustSizes<-rep(0,chosenNClusters)
      for(c in 1:chosenNClusters){
         clustSizes[c]<-length(which(clustVec==c))
      }

      clusteringPred<-NULL
      if(nPredictSubjects>0){
         clusteringPred<-match(optZPredict,uZFit)
      }
      avgSilhouetteWidth<-NULL
      close(zFile)
      
	}else{
   
		if(is.null(maxNClusters)){
			# Determine the maximum number of clusters
			nMembersFileName<-file.path(directoryPath,paste(fileStem,'_nMembers.txt',sep=''))
			nMembersFile<-file(nMembersFileName,open="r")
			nClustersFileName<- file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
			nClustersFile<-file(nClustersFileName,open="r")
   	
			# Restrict to sweeps after burn in
			firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
			lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter
			maxNClusters<-0
   	
			for(sweep in firstLine:lastLine){
				if(sweep==firstLine){
					skipVal<-firstLine-1
				}else{
					skipVal<-0
				}
	   
				# Get the current number of members for each cluster
				nClusters<-scan(nClustersFile,what=integer(),skip=skipVal,n=1,quiet=T)	
				currNMembers<-scan(nMembersFile,what=integer(),skip=skipVal,n=nClusters+1,quiet=T)
				currNMembers<-currNMembers[1:nClusters]
				# Find the number of non-empty clusters
				nNotEmpty<-sum(currNMembers>0)
				if(nNotEmpty>maxNClusters){
					maxNClusters<-nNotEmpty
				}
			}   
			# Add on another 5 just to make sure bound is safe (but don't let it exceed no. of subjects -1)
			maxNClusters<-min(maxNClusters+5,nSubjects-1)

			close(nMembersFile)
			close(nClustersFile)
		}   
   
		# If the input was a list of dissimilarity matrices then take the average
		if(is.list(disSimMat)){
			for(i in 1:length(disSimMat)){
				if(i==1){
					tmpMat<-disSimMat[[i]]      
				}else{
					tmpMat<-tmpMat+disSimMat[[i]]
				}
			}
			tmpMat<-tmpMat/length(disSimMat)
			disSimMat<-tmpMat
		}
   
		# Loop over the possible number of clusters
		avgSilhouetteWidth<--1.0;
		cat(paste("Max no of possible clusters:",maxNClusters,"\n"))
		for(c in 2:maxNClusters){
			cat(paste("Trying",c,"clusters\n"))
			tmpObj<-pam(disSimMat,k=c,diss=T)
			# Check whether the silhouette width from this clustering improves previous best
			if(avgSilhouetteWidth<tmpObj$silinfo$avg.width){
				avgSilhouetteWidth<-tmpObj$silinfo$avg.width
				chosenNClusters<-c
				clustVec<-tmpObj$clustering
				clustSizes<-tmpObj$clusinfo[,1]
				# The id of the objects chosen as the medoids
				clustMedoids<-tmpObj$id.med
			}
		}
   
		# Work out the clustering of the prediction objects
		clusteringPred<-NULL
		if(nPredictSubjects>0){
			disSimMatPred<-matrix(disSimMatPred,nrow=nPredictSubjects,byrow=T)
			clusteringPred<-rep(0,nPredictSubjects)
			for(i in 1:nPredictSubjects){
				tmpVec<-disSimMatPred[i,clustMedoids]
				whichMin <- which(tmpVec==min(tmpVec))
				if (length(whichMin)>1) {
					clusteringPred[i]<-sample(whichMin,1)
				} else {
					clusteringPred[i]<-whichMin
				}
			}
		}
	}
	
	return(list("clusObjRunInfoObj"=disSimObj$disSimRunInfoObj,
		"clusObjDisSimMat"=disSimObj$disSimMat,
		"nClusters"=chosenNClusters,
		"clusterSizes"=clustSizes,
		"clustering"=clustVec,
		"avgSilhouetteWidth"=avgSilhouetteWidth,
		"clusteringPred"=clusteringPred))
}

# Function to take the optimal clustering and computing the risk and probability
# profile
calcAvgRiskAndProfile<-function(clusObj,includeFixedEffects=F,proportionalHazards=F){

	clusObjRunInfoObj=NULL
	directoryPath=NULL
	fileStem=NULL
	xModel=NULL
	nCategories=NULL
	varSelect=NULL
	varSelectType=NULL
	includeResponse=NULL
	nFixedEffects=NULL
	reportBurnIn=NULL
	nBurn=NULL
	nFilter=NULL
	nSweeps=NULL
	nClusters=NULL
	clustering=NULL
	nCategoriesY=NULL
	nCovariates=NULL
	nDiscreteCovs=NULL
	nContinuousCovs=NULL
	nSubjects=NULL
	nPredictSubjects=NULL
	yModel=NULL
	wMat=NULL
	yMat=NULL
	weibullFixedShape=NULL

	for (i in 1:length(clusObj)) assign(names(clusObj)[i],clusObj[[i]])
	for (i in 1:length(clusObjRunInfoObj)) assign(names(clusObjRunInfoObj)[i],clusObjRunInfoObj[[i]])
	 
	# Construct the number of clusters file name
	nClustersFileName<-file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
	nClustersFile<-file(nClustersFileName,open="r")
	
	# Construct the allocation file name
	zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
	zFile<-file(zFileName,open="r")
	
	if(xModel=="Discrete"){
		# Construct the allocation file name
		phiFileName <- file.path(directoryPath,paste(fileStem,'_phi.txt',sep=''))
		phiFile<-file(phiFileName,open="r")
		# Get the maximum number of categories
		maxNCategories<-max(nCategories)
		if(varSelect){
			nullPhiFileName <- file.path(directoryPath,paste(fileStem,'_nullPhi.txt',sep=''))
			if(varSelectType=="Continuous"){
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}else{
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_gamma.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}
		}
	
	}else if(xModel=="Normal"){
		muFileName <- file.path(directoryPath,paste(fileStem,'_mu.txt',sep=''))
		muFile<-file(muFileName,open="r")
		SigmaFileName <- file.path(directoryPath,paste(fileStem,'_Sigma.txt',sep=''))
		SigmaFile<-file(SigmaFileName,open="r")
		if(varSelect){
			nullMuFileName <- file.path(directoryPath,paste(fileStem,'_nullMu.txt',sep=''))
			if(varSelectType=="Continuous"){
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}else{
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_gamma.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}
		}
	}else if(xModel=="Mixed"){
		# Construct the allocation file name
		phiFileName <- file.path(directoryPath,paste(fileStem,'_phi.txt',sep=''))
		phiFile<-file(phiFileName,open="r")
		# Get the maximum number of categories
		maxNCategories<-max(nCategories)
		muFileName <- file.path(directoryPath,paste(fileStem,'_mu.txt',sep=''))
		muFile<-file(muFileName,open="r")
		SigmaFileName <- file.path(directoryPath,paste(fileStem,'_Sigma.txt',sep=''))
		SigmaFile<-file(SigmaFileName,open="r")
		if(varSelect){
			nullPhiFileName <- file.path(directoryPath,paste(fileStem,'_nullPhi.txt',sep=''))
			nullMuFileName <- file.path(directoryPath,paste(fileStem,'_nullMu.txt',sep=''))
			if(varSelectType=="Continuous"){
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}else{
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_gamma.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}
		}
	}
	
	if(includeResponse){
		# Construct the theta file name
		thetaFileName <- file.path(directoryPath,paste(fileStem,'_theta.txt',sep=''))
		thetaFile<-file(thetaFileName,open="r")
		if (yModel=="Survival"&&!weibullFixedShape){
			nuFileName <- file.path(directoryPath,paste(fileStem,'_nu.txt',sep=''))
			nuFile<-file(nuFileName,open="r")
		}
		if(nFixedEffects>0){
			# Construct the fixed effect coefficient file name
			betaFileName <-file.path(directoryPath,paste(fileStem,'_beta.txt',sep=''))
			betaFile<-file(betaFileName,open="r")
		}
		if (yModel=="Survival"&&weibullFixedShape){
				nuFileName<-file.path(directoryPath,paste(fileStem,'_nu.txt',sep=''))
				nuFile<-file(nuFileName,open="r")
				nu<-(read.table(nuFile)[,1])
				close(nuFile)
		}
	} 
	
	# Restrict to sweeps after burn in
	firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
	lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter
	nSamples<-lastLine-firstLine+1
	
	# Make a list of the subjects in each of the optimal clusters
	optAlloc<-vector("list",nClusters)
	for(c in 1:nClusters){
		optAlloc[[c]]<-which(clustering==c)
	}
	
	if(includeResponse){
		# Initialise the object for storing the risks
		riskArray<-array(0,dim=c(nSamples,nClusters,nCategoriesY))
		thetaArray<-array(0,dim=c(nSamples,nClusters,nCategoriesY))
		if (yModel=="Survival"&&!weibullFixedShape) nuArray<-array(0,dim=c(nSamples,nClusters))
		if(nFixedEffects>0){
			betaArray<-array(0,dim=c(nSamples,nFixedEffects,nCategoriesY))
		}
	}else{
		riskArray<-NULL
	}

	# Initialise the object for storing the profiles
	if(xModel=='Discrete'){
		phiArray<-array(dim=c(nSamples,nClusters,nCovariates,maxNCategories))
		if(varSelect){
			phiStarArray<-array(dim=c(nSamples,nClusters,nCovariates,maxNCategories))
			tmpCurrNullPhi<-scan(nullPhiFileName,what=double(),quiet=T)
			tmpCurrNullPhi<-array(tmpCurrNullPhi,dim=c(maxNCategories,nCovariates))
			currNullPhi<-array(dim=c(1,maxNCategories,nCovariates))
			currNullPhi[1,,]<-tmpCurrNullPhi
		}else{
			phiStarArray<-NULL
		}
	}else if(xModel=='Normal'){
		muArray<-array(dim=c(nSamples,nClusters,nCovariates))
		if(varSelect){
			muStarArray<-array(dim=c(nSamples,nClusters,nCovariates))
			currNullMu<-scan(nullMuFileName,what=double(),quiet=T)
			currNullMu<-array(currNullMu,dim=c(nCovariates,1))
			currNullMu<-t(currNullMu)
		}else{
			muStarArray<-NULL
		}
		sigmaArray<-array(dim=c(nSamples,nClusters,nCovariates,nCovariates))
	}else if(xModel=='Mixed'){
		phiArray<-array(dim=c(nSamples,nClusters,nDiscreteCovs,maxNCategories))
		if(varSelect){
			phiStarArray<-array(dim=c(nSamples,nClusters,nDiscreteCovs,maxNCategories))
			tmpCurrNullPhi<-scan(nullPhiFileName,what=double(),quiet=T)
			tmpCurrNullPhi<-array(tmpCurrNullPhi,dim=c(maxNCategories,nDiscreteCovs))
			currNullPhi<-array(dim=c(1,maxNCategories,nDiscreteCovs))
			currNullPhi[1,,]<-tmpCurrNullPhi
		}else{
			phiStarArray<-NULL
		}
		muArray<-array(dim=c(nSamples,nClusters,nContinuousCovs))
		if(varSelect){
			muStarArray<-array(dim=c(nSamples,nClusters,nContinuousCovs))
			currNullMu<-scan(nullMuFileName,what=double(),quiet=T)
			currNullMu<-array(currNullMu,dim=c(nContinuousCovs,1))
			currNullMu<-t(currNullMu)
		}else{
			muStarArray<-NULL
		}
		sigmaArray<-array(dim=c(nSamples,nClusters,nContinuousCovs,nContinuousCovs))
	}
	

	for(sweep in firstLine:lastLine){
		if(sweep==firstLine){
			skipVal<-firstLine-1
		}else{
			skipVal<-0
		}
	
		if(sweep-firstLine==0||(sweep-firstLine+1)%%1000==0){
			cat(paste("Processing sweep",sweep-firstLine+1,"of ",lastLine-firstLine+1,"\n"))
		}
		currMaxNClusters<-scan(nClustersFile,what=integer(),skip=skipVal,n=1,quiet=T)
	
		# Get the current allocation data for this sweep
		currZ<-scan(zFile,what=integer(),skip=skipVal,n=nSubjects+nPredictSubjects,quiet=T)
		currZ<-1+currZ

		if(includeResponse){
			# Get the risk data corresponding to this sweep
			if (yModel=="Categorical") {
				currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currMaxNClusters*(nCategoriesY-1),quiet=T)
				currTheta<-matrix(currThetaVector,ncol=(nCategoriesY-1),byrow=T)
				currTheta<-cbind(rep(0,dim(currTheta)[1]),currTheta)
			} else {
				currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCategoriesY,quiet=T)
				currTheta<-matrix(currThetaVector,ncol=nCategoriesY,byrow=T)
				if (yModel=="Survival"&&!weibullFixedShape){
					currNuVector<-scan(nuFile,what=double(),skip=skipVal,n=currMaxNClusters,quiet=T)
					currNu<-c(currNuVector)
				}
			}
			if(nFixedEffects>0){
				if (yModel=="Categorical") {
					currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects*(nCategoriesY-1),quiet=T)
					currBeta<-matrix(currBetaVector,ncol=(nCategoriesY-1),byrow=T)
					currBeta<-cbind(rep(0,dim(currBeta)[1]),currBeta)
				} else {
					currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects*nCategoriesY,quiet=T)
					currBeta<-matrix(currBetaVector,ncol=nCategoriesY,byrow=T)
				}
				betaArray[sweep-firstLine+1,,]<-currBeta
			}
			# Calculate the average risk (over subjects) for each cluster

			for(c in 1:nClusters){
				currLambdaVector<-currTheta[currZ[optAlloc[[c]]],]
				currLambda<-matrix(currLambdaVector,ncol=nCategoriesY)
				if(includeFixedEffects&&nFixedEffects>0){
					if (yModel=="Categorical"){
						for (i in 1:length(optAlloc[[c]])){
							for (k in 1:nCategoriesY) currLambda[i,k]<-currLambda[i,k]+
								t(as.matrix(wMat[optAlloc[[c]][i],]))%*%
								currBeta[,yMat[optAlloc[[c]][i]]+1]
						}	
					} else {
						currLambda<-currLambda+as.matrix(wMat[optAlloc[[c]],])%*%currBeta
					}
				}
				if(yModel=="Poisson"){
					currRisk<-exp(currLambda )
				}else if(yModel=="Bernoulli"||yModel=="Binomial"){
					currRisk<-1.0/(1.0+exp(-currLambda))
				}else if(yModel=="Normal"){
					currRisk<-currLambda
				}else if(yModel=="Categorical"){
					currRisk<-matrix(0,ncol=length(optAlloc[[c]]),nrow=nCategoriesY)
					currRisk<-exp(currLambda)/rowSums(exp(currLambda))
				}else if(yModel=="Survival"){
					if (proportionalHazards){
						currRisk<-exp(currLambda)
					}else{
						if (!weibullFixedShape){
							currNuVector<-currNu[currZ[optAlloc[[c]]]]
						} else {
							currNuVector<-nu[sweep]
						}
						currRisk<-1/((exp(currLambda))^(1/currNuVector))*gamma(1+1/currNuVector)
					}

				}
				riskArray[sweep-firstLine+1,c,]<-apply(currRisk,2,mean)
				thetaArray[sweep-firstLine+1,c,]<-apply(as.matrix(currTheta[currZ[optAlloc[[c]]],],ncol=nCategoriesY),2,mean)
				if(yModel=="Survival"&&!weibullFixedShape){
					nuArray[sweep-firstLine+1,c]<-mean(currNuVector)
				}
			}
		}

		# Calculate the average profile (over subjects) for each cluster
		if(xModel=='Discrete'){
			currPhi<-scan(phiFile,what=double(),
				skip=skipVal,n=currMaxNClusters*maxNCategories*nCovariates,quiet=T)
			# This is slightly convoluted, because of the way that R reads in by column
			# I switched the order of categories and covariates in column below, and then
			# take the transpose to correct in the loop
			currPhi<-array(currPhi,dim=c(currMaxNClusters,maxNCategories,nCovariates))
			if(varSelect){
				# We increase dimensions of currGamma using duplicates, to
				# enable easier calculation of phiStar
				if(varSelectType=='BinaryCluster'){
					tmpCurrGamma<-scan(gammaFile,what=double(),
						skip=skipVal,n=nCovariates*currMaxNClusters,quiet=T)
					tmpCurrGamma<-array(tmpCurrGamma,dim=c(currMaxNClusters,nCovariates))
				}else{
					tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=nCovariates,quiet=T)
					tmpCurrGamma<-array(tmpCurrGamma,dim=c(nCovariates,currMaxNClusters))
					tmpCurrGamma<-t(tmpCurrGamma)
				}
				currGamma<-array(dim=c(currMaxNClusters,maxNCategories,nCovariates))
				for(p in 1:maxNCategories){
					currGamma[,p,]<-tmpCurrGamma
				}
			}	
			for(c in 1:nClusters){
				phiArray[sweep-firstLine+1,c,,]<-t(apply(array(currPhi[currZ[optAlloc[[c]]],,],
					dim=c(length(optAlloc[[c]]),
					dim(currPhi)[2],dim(currPhi)[3])),2:3,mean))
				if(varSelect){
					phiStarArray[sweep-firstLine+1,c,,]<-t(apply(array(currGamma[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currGamma)[2],
						dim(currGamma)[3]))*array(currPhi[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currPhi)[2],dim(currPhi)[3]))+
						(1-array(currGamma[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currGamma)[2],dim(currGamma)[3])))*
						array(currNullPhi[rep(1,length(optAlloc[[c]])),,],
						dim=c(length(optAlloc[[c]]),dim(currNullPhi)[2],
						dim(currNullPhi)[3])),2:3,mean))

				}
			}
		}else if(xModel=='Normal'){
			# mu stored like phi
			currMu<-scan(muFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates,quiet=T)
			currMu<-array(currMu,dim=c(currMaxNClusters,nCovariates))
			if(varSelect){
				# We increase dimensions of nullPhi and currGamma using duplicates, to
				# enable easier calculation of phiStar
				# We increase dimensions of currGamma using duplicates, to
				# enable easier calculation of phiStar
				if(varSelectType=='BinaryCluster'){
					tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates,quiet=T)
					currGamma<-array(tmpCurrGamma,dim=c(currMaxNClusters,nCovariates))
				}else{
					tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=nCovariates,quiet=T)
					tmpCurrGamma<-array(tmpCurrGamma,dim=c(nCovariates,currMaxNClusters))
					currGamma<-t(tmpCurrGamma)
				}
	
			}
			for(c in 1:nClusters){
				muArray[sweep-firstLine+1,c,]<-apply(matrix(currMu[currZ[optAlloc[[c]]],],ncol=nCovariates),2,mean)
				if(varSelect){
					muStarArray[sweep-firstLine+1,c,]<-apply(matrix(currGamma[currZ[optAlloc[[c]]],],
						ncol=nCovariates)*matrix(currMu[currZ[optAlloc[[c]]],],ncol=nCovariates)+
						matrix(1-currGamma[currZ[optAlloc[[c]]],],ncol=nCovariates)*
						matrix(currNullMu[rep(1,length(optAlloc[[c]])),],ncol=nCovariates),2,mean)
				}
			}
	
			currSigma<-scan(SigmaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates*nCovariates,quiet=T)
			currSigma<-array(currSigma,dim=c(currMaxNClusters,nCovariates,nCovariates))
			for(c in 1:nClusters){
				sigmaArray[sweep-firstLine+1,c,,]<-apply(array(currSigma[currZ[optAlloc[[c]]],,],
					dim=c(length(optAlloc[[c]]),dim(currSigma)[2],dim(currSigma)[3])),2:3,mean)
			}
		}else if(xModel=='Mixed'){
			currPhi<-scan(phiFile,what=double(),
				skip=skipVal,n=currMaxNClusters*maxNCategories*nDiscreteCovs,quiet=T)
			# This is slightly convoluted, because of the way that R reads in by column
			# I switched the order of categories and covariates in column below, and then
			# take the transpose to correct in the loop
			currPhi<-array(currPhi,dim=c(currMaxNClusters,maxNCategories,nDiscreteCovs))
			# mu stored like phi
			currMu<-scan(muFile,what=double(),skip=skipVal,n=currMaxNClusters*nContinuousCovs,quiet=T)
			currMu<-array(currMu,dim=c(currMaxNClusters,nContinuousCovs))
			if(varSelect){
				# We increase dimensions of currGamma using duplicates, to
				# enable easier calculation of phiStar
				if(varSelectType=='BinaryCluster'){
					tmpCurrGamma<-scan(gammaFile,what=double(),
						skip=skipVal,n=nCovariates*currMaxNClusters,quiet=T)
					tmpCurrGamma<-array(tmpCurrGamma,dim=c(currMaxNClusters,nCovariates))
				}else{
					tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=nCovariates,quiet=T)
					tmpCurrGamma<-array(tmpCurrGamma,dim=c(nCovariates,currMaxNClusters))
					tmpCurrGamma<-t(tmpCurrGamma)
				}
				currGamma<-array(dim=c(currMaxNClusters,maxNCategories,nCovariates))
				for(p in 1:maxNCategories){
					currGamma[,p,]<-tmpCurrGamma
				}
			}	
			for(c in 1:nClusters){
				phiArray[sweep-firstLine+1,c,,]<-t(apply(array(currPhi[currZ[optAlloc[[c]]],,],
					dim=c(length(optAlloc[[c]]),
					dim(currPhi)[2],dim(currPhi)[3])),2:3,mean))
				if(varSelect){
					currGammaD<-array(currGamma[,,1:nDiscreteCovs],dim=c(currMaxNClusters,maxNCategories,nDiscreteCovs))
					phiStarArray[sweep-firstLine+1,c,,]<-t(apply(array(currGammaD[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currGammaD)[2],
						dim(currGammaD)[3]))*array(currPhi[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currPhi)[2],dim(currPhi)[3]))+
						(1-array(currGammaD[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currGammaD)[2],dim(currGammaD)[3])))*
						array(currNullPhi[rep(1,length(optAlloc[[c]])),,],
						dim=c(length(optAlloc[[c]]),dim(currNullPhi)[2],
						dim(currNullPhi)[3])),2:3,mean))
				}
				muArray[sweep-firstLine+1,c,]<-apply(matrix(currMu[currZ[optAlloc[[c]]],],ncol=nContinuousCovs),2,mean)
				if(varSelect){
					currGammaC<-array(currGamma[,1,(nDiscreteCovs+1):nCovariates],
						dim=c(currMaxNClusters,nContinuousCovs))
					muStarArray[sweep-firstLine+1,c,]<-apply(matrix(currGammaC[currZ[optAlloc[[c]]],],
						ncol=nContinuousCovs)*matrix(currMu[currZ[optAlloc[[c]]],],ncol=nContinuousCovs)+
						matrix(1-currGammaC[currZ[optAlloc[[c]]],],ncol=nContinuousCovs)*
						matrix(currNullMu[rep(1,length(optAlloc[[c]])),],ncol=nContinuousCovs),2,mean)
				}
			}
	
			currSigma<-scan(SigmaFile,what=double(),skip=skipVal,n=currMaxNClusters*nContinuousCovs*nContinuousCovs,quiet=T)
			currSigma<-array(currSigma,dim=c(currMaxNClusters,nContinuousCovs,nContinuousCovs))
			for(c in 1:nClusters){
				sigmaArray[sweep-firstLine+1,c,,]<-apply(array(currSigma[currZ[optAlloc[[c]]],,],
					dim=c(length(optAlloc[[c]]),dim(currSigma)[2],dim(currSigma)[3])),2:3,mean)
			}
		}
	}

	# Calculate the empiricals
	empiricals<-rep(0,nClusters)
	if(!is.null(yModel)){
		for(c in 1:nClusters){
			if(yModel=='Bernoulli'||yModel=='Normal'||yModel=='Survival'){
				empiricals[c]<-mean(yMat[optAlloc[[c]],1])
			}else if(yModel=='Binomial'){
				empiricals[c]<-mean(yMat[optAlloc[[c]],1]/yMat[optAlloc[[c]],2])
			}else if(yModel=='Poisson'){
				empiricals[c]<-mean(yMat[optAlloc[[c]],1]/yMat[optAlloc[[c]],2])
			#}else if(yModel=='Categorical'){
			# no empiricals for categorical outcome				
			}
		}
	}

	if(xModel=='Discrete'){
		out<-list('riskProfClusObj'=clusObj,'risk'=riskArray,'profile'=phiArray,'profileStar'=phiStarArray,'empiricals'=empiricals)
	}else if(xModel=='Normal'){
		out<-list('riskProfClusObj'=clusObj,'risk'=riskArray,
			'profile'=muArray,'profileStar'=muStarArray,
			'profileStdDev'=sigmaArray,'empiricals'=empiricals)
	}else if(xModel=='Mixed'){
		out<-list('riskProfClusObj'=clusObj,'risk'=riskArray,
			'profilePhi'=phiArray,'profileStarPhi'=phiStarArray,
			'profileMu'=muArray,'profileStarMu'=muStarArray,
			'profileStdDev'=sigmaArray,'empiricals'=empiricals)
	}
	
	close(zFile)
	close(nClustersFile)
	if(xModel=="Discrete"){
		close(phiFile)
		if(varSelect){
			close(gammaFile)
		}
	}else if(xModel=="Normal"){
		close(muFile)
		close(SigmaFile)
		if(varSelect){
			close(gammaFile)
		}
	}else if(xModel=="Mixed"){
		close(phiFile)
		close(muFile)
		close(SigmaFile)
		if(varSelect){
			close(gammaFile)
		}
	}

	
	if(includeResponse){
		close(thetaFile)
		if(nFixedEffects>0){
			close(betaFile)
		}
		if (yModel=="Survival"&&!weibullFixedShape) {
			out$nuArray<-nuArray
			close(nuFile)
		}
	}
	
	return(out)
}
	
		

# Plot output values
plotRiskProfile<-function(riskProfObj,outFile,showRelativeRisk=F,orderBy=NULL,whichClusters=NULL,whichCovariates=NULL,useProfileStar=F,riskLim=NULL){

	riskProfClusObj=NULL
	clusObjRunInfoObj=NULL
	includeResponse=NULL
	yModel=NULL
	profileStar=NULL
	xModel=NULL
	whicCov=NULL
	nCategoriesY=NULL
	cluster=NULL
	prob=NULL
	meanProb=NULL
	fillColor=NULL
	lowerProb=NULL
	upperProb=NULL
	meanRisk=NULL
	lowerRisk=NULL
	upperRisk=NULL
	clusterSize=NULL
	mu=NULL
	meanMu=NULL
	lowerMu=NULL
	upperMu=NULL
	sigma=NULL
	meanSigma=NULL
	lowerSigma=NULL
	upperSigma=NULL
	weibullFixedShape=NULL
	nu=NULL
	meanNu=NULL
	lowerNu=NULL
	upperNu=NULL	

	for (i in 1:length(riskProfObj)) assign(names(riskProfObj)[i],riskProfObj[[i]])
	for (i in 1:length(riskProfClusObj)) assign(names(riskProfClusObj)[i],riskProfClusObj[[i]])
	for (i in 1:length(clusObjRunInfoObj)) assign(names(clusObjRunInfoObj)[i],clusObjRunInfoObj[[i]])

	if (nClusters==1) stop("Cannot produce plots because only one cluster has been found.")

	if(includeResponse){
		if(yModel=="Normal"){
			showRelativeRisk<-F
		}
	}

	if(useProfileStar){
		profile<-profileStar
	}
	if(!is.null(whichCovariates)){
		if (!is.numeric(whichCovariates)){
			whichCovariatesTmp<-vector()
			for (k in 1:length(whichCovariates)){
				whichCovariatesTmp[k]<-which(riskProfClusObj$clusObjRunInfoObj$covNames==whichCovariates[k])
			}
			whichCovariates<-whichCovariatesTmp
		}
		if(xModel=='Discrete'){
			profile<-profile[,,whichCovariates,]
			nCategories<-nCategories[whichCovariates]
			covNames<-covNames[whichCovariates]
			nCovariates<-length(whichCovariates)
		}else if(xModel=='Normal'){
			profile<-profile[,,whichCovariates]
			profileStdDev<-profileStdDev[,,whichCovariates,whichCovariates]
			covNames<-covNames[whichCovariates]
			nCovariates<-length(whichCovariates)
		}else if(xModel=='Mixed'){
			nDiscreteCovsAll <- nDiscreteCovs
			nContinuousCovsAll <- nContinuousCovs
			whichDiscreteCovs <- whichCovariates[which(whichCovariates<=nDiscreteCovs)]
			whichContinuousCovs <- whichCovariates[which(whichCovariates>nDiscreteCovs)]
			discreteCovs <- discreteCovs[whichDiscreteCovs]
			nDiscreteCovs <- length(discreteCovs)
			continuousCovs <- continuousCovs[whichContinuousCovs-nDiscreteCovsAll]
			nContinuousCovs <- length(continuousCovs)
			profilePhi<-profilePhi[,,whichDiscreteCovs,]
			nCategories<-nCategories[whichDiscreteCovs]
			profileMu<-profileMu[,,whichContinuousCovs-nDiscreteCovsAll]
			profileStdDev<-profileStdDev[,,whichContinuousCovs-nDiscreteCovsAll,whichContinuousCovs-nDiscreteCovsAll]
			covNames<-c(discreteCovs,continuousCovs)
			nCovariates<-length(covNames)

		}
	}

	png(outFile,width=1200,height=800)
	orderProvided<-F
	
	if(!is.null(orderBy)){
		if(!includeResponse){
			if(orderBy!='Empirical'&&orderBy!='ClusterSize'&&!orderBy%in%covNames){
				if(is.numeric(orderBy)){
					if(length(orderBy)==nClusters){
						orderProvided<-T
						meanSortIndex<-orderBy
					}else{
						cat("Order vector provided not of same length as number of clusters. Reverting to default ordering.\n")	
						orderBy<-NULL
					}
				orderBy<-NULL
				}
				#orderBy<-NULL
			}
		}else{
			if(orderBy!='Risk'&&orderBy!='Empirical'&&orderBy!='ClusterSize'&&!orderBy%in%covNames){
				if(is.numeric(orderBy)){
					if(length(orderBy)==nClusters){
						orderProvided<-T
						meanSortIndex<-orderBy
					}else{
						cat("Order vector provided not of same length as number of clusters. Reverting to default ordering.\n")	
						orderBy<-NULL
					}
				orderBy<-NULL
				}
				#orderBy<-NULL
			}
		}

	}
	
	# Set up the layout for the plot
	plotLayout<-grid.layout(ncol = nCovariates+2, nrow = 6)
	grid.newpage()
	pushViewport(viewport(layout = plotLayout))

	if(!orderProvided){
		if(!is.null(risk)){
			if(is.null(orderBy)){
				# Default is to order by posterior theta risk
				# Compute the means
				orderStat<-apply(risk,2,median)
			}else{
				if(orderBy=='Risk'){
					orderStat<-apply(risk,2,median)
				}else if(orderBy=='Empirical'){
					orderStat<-empiricals
				}else if(orderBy=='ClusterSize'){
					orderStat<-clusterSizes
				}else{
					whichCov<-match(orderBy,covNames)
					if(xModel=='Normal'){
						orderStat<-apply(profile[,,whichCov],2,median)
					}else{
						# This assumes that there is some order to the categories
						# and then uses an expected value
						tmpMat<-profile[,,whichCov,1]
						if(nCategories[whichCov]>1){
							for(k in 2:nCategories[whichCov]){
								tmpMat<-tmpMat+k*profile[,,whichCov,k]
							}
						}
						orderStat<-apply(tmpMat,2,median)
					}
				}
			}
		}else{
			if(is.null(orderBy)){
				# Default is to order by empirical risk
				orderStat<-empiricals
			}else{
				if(orderBy=='Empirical'){
					orderStat<-empiricals
				}else if(orderBy=='ClusterSize'){
					orderStat<-clusterSizes
				}else{
					whichCov<-match(orderBy,covNames)
					if(xModel=='Normal'){
						orderStat<-apply(profile[,,whichCov],2,median)
					}else{
						# This assumes that there is some order to the categories
						# and then uses an expected value
						tmpMat<-profile[,,whichCov,1]
						if(nCategories[whichCov]>1){
							for(k in 2:(nCategories[whichCov])){
								tmpMat<-tmpMat+k*profile[,,whichCov,k]
							}
						}
						orderStat<-apply(tmpMat,2,median)
					}
				}
			}	
		}
		# Sort into ascending mean size
		meanSortIndex<-order(orderStat,decreasing=F)
	}
	if(includeResponse){
		# Reorder the risk matrix
		riskDim<-dim(risk)
		risk<-array(risk[,meanSortIndex,],dim=riskDim)
		if(showRelativeRisk){
			for(c in nClusters:1){
				risk[,c,]<-risk[,c,]/risk[,1,]
			}
		}
		# reorder the nu matrix 
		if (yModel=="Survival"&&!weibullFixedShape){
			nuDim<-dim(nuArray)
			nuArray<-array(nuArray[,meanSortIndex],dim=nuDim)
		}
	}

	# Reorder the cluster sizes
	clusterSizes<-clusterSizes[meanSortIndex]
	# Reorder the empiricals
	empiricals<-empiricals[meanSortIndex]
	meanEmpirical<-sum(empiricals*clusterSizes)/sum(clusterSizes)

	if(includeResponse){
		# Recompute the means and now also credible intervals
		riskMeans<-apply(risk,2,mean,trim=0.005)
		riskMean<-sum(riskMeans*clusterSizes)/sum(clusterSizes)
		riskLower<-apply(risk,2,quantile,0.05)
		riskUpper<-apply(risk,2,quantile,0.95)
		# The next line is to avoid outliers spoiling plot scales
		plotMax<-max(riskUpper)
		
		# Get the plot colors
		riskColor<-ifelse(riskLower>rep(riskMean,nClusters),"high",
		ifelse(riskUpper<rep(riskMean,nClusters),"low","avg"))
		if (yModel=="Categorical"){		
			riskDF<-data.frame("risk"=c(),"category"=c(),"cluster"=c(),"meanRisk"=c(),
				"lowerRisk"=c(),"upperRisk"=c(),"fillColor"=c())
		} else {
			riskDF<-data.frame("risk"=c(),"cluster"=c(),"meanRisk"=c(),
				"lowerRisk"=c(),"upperRisk"=c(),"fillColor"=c())
		}
		if (yModel=="Survival"&&!weibullFixedShape){
			# Recompute the means and now also credible intervals
			nuMeans<-apply(nuArray,2,mean,trim=0.005)
			nuMean<-sum(nuMeans*clusterSizes)/sum(clusterSizes)
			nuLower<-apply(nuArray,2,quantile,0.05)
			nuUpper<-apply(nuArray,2,quantile,0.95)

			nuDF<-data.frame("nu"=c(),"cluster"=c(),"meanNu"=c(),
				"lowerNu"=c(),"upperNu"=c(),"fillColor"=c())
		}

	}else{
		riskColor<-ifelse(empiricals>rep(meanEmpirical,length(empiricals)),"high",
		ifelse(empiricals<rep(meanEmpirical,nClusters),"low","avg"))
	}


	if(is.null(whichClusters)){
		whichClusters<-1:nClusters
	}
	nClusters<-length(whichClusters)
	
	empiricalDF<-data.frame("empiricals"=c(),"meanEmpirical"=c(),"cluster"=c(),"fillColor"=c())
	sizeDF<-data.frame("clusterSize"=c(),"cluster"=c(),"fillColor"=c())
	
	# Restructure the data for plotting
	for(c in whichClusters){
		if(includeResponse){
			if (yModel=="Categorical"){
				plotRisk<-risk[,c,]
				nPoints<-dim(plotRisk)[1]
				for (k in 1:nCategoriesY){				
					riskDF<-rbind(riskDF,data.frame("risk"=plotRisk[,k],
						"category"=rep(k,nPoints),
						"cluster"=rep(c,nPoints),
						"meanRisk"=rep(riskMean,nPoints),
						"lowerRisk"=rep(riskLower[c],nPoints),
						"upperRisk"=rep(riskUpper[c],nPoints),
						"fillColor"=rep(riskColor[c],nPoints)))
				}
			} else {
				plotRisk<-risk[,c,]
				plotRisk<-plotRisk[plotRisk<plotMax]
				nPoints<-length(plotRisk)
				riskDF<-rbind(riskDF,data.frame("risk"=plotRisk,"cluster"=rep(c,nPoints),
					"meanRisk"=rep(riskMean,nPoints),
					"lowerRisk"=rep(riskLower[c],nPoints),
					"upperRisk"=rep(riskUpper[c],nPoints),
					"fillColor"=rep(riskColor[c],nPoints)))
			}
			if (yModel=="Survival"&&!weibullFixedShape){
				plotNu<-nuArray[,c]
				nPoints<-length(plotNu)
				nuDF<-rbind(nuDF,data.frame("nu"=plotNu,"cluster"=rep(c,nPoints),
					"meanNu"=rep(nuMean,nPoints),
					"lowerNu"=rep(nuLower[c],nPoints),
					"upperNu"=rep(nuUpper[c],nPoints),
					"fillColor"=rep(riskColor[c],nPoints)))
			}
		}
		empiricalDF<-rbind(empiricalDF,
			data.frame("empiricals"=empiricals[c],
			"meanEmpirical"=meanEmpirical,"cluster"=c,"fillColor"=riskColor[c]))
		sizeDF<-rbind(sizeDF,
			data.frame("clusterSize"=clusterSizes[c],"cluster"=c,"fillColor"=riskColor[c]))
	}

	if(includeResponse){
		if(yModel=='Categorical'){
			riskDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
				"lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
			for(k in 1:nCategoriesY){
		
				probMat<-risk[,,k]
				nPoints<-nrow(probMat)
				probMeans<-apply(probMat,2,mean)
				probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
				probLower<-apply(probMat,2,quantile,0.05)
				probUpper<-apply(probMat,2,quantile,0.95)
		
				# Get the plot colors
				probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
				ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))

				for(c in whichClusters){
					riskDF<-rbind(riskDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
						"category"=rep(k-1,nPoints),
						"meanProb"=rep(probMean,nPoints),
						"lowerProb"=rep(probLower[c],nPoints),
						"upperProb"=rep(probUpper[c],nPoints),
						"fillColor"=rep(probColor[c],nPoints)))
						 rownames(riskDF)<-seq(1,nrow(riskDF),1)
				
				}
			}

			plotObj<-ggplot(riskDF)
			plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
			plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
			# Margin order is (top,right,bottom,left)
			plotObj<-plotObj+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))+
				theme(plot.margin=unit(c(0,0,0,0),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=2))
#		}else if (yModel=="Survival"&&!weibullFixedShape){
#			rownames(riskDF)<-seq(1,nrow(riskDF),by=1)
#
#			# Create the risk plot
#			plotObj<-ggplot(riskDF)
#			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=risk,yintercept=meanRisk))
#			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=risk,fill=as.factor(fillColor)),outlier.size=0.5)
#			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerRisk,colour=as.factor(fillColor)),size=1.5)
#			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperRisk,colour=as.factor(fillColor)),size=1.5)
#			plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
#				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
#				theme(legend.position="none")+
#				labs(x="Cluster",y=ifelse(showRelativeRisk,'RR',
#				ifelse(yModel=="Categorical"||yModel=="Bernoulli"||yModel=="Binomial","Probability","E[Y]")))
#			plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
#			plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
#			# Margin order is (top,right,bottom,left)
#			plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
#			print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=2))
#
#			rownames(nuDF)<-seq(1,nrow(nuDF),by=1)
#			# Create the nu plot
#			plotObj<-ggplot(nuDF)
#			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=nu,yintercept=meanNu))
#			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=nu,fill=as.factor(fillColor)),outlier.size=0.5)
#			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerNu,colour=as.factor(fillColor)),size=1.5)
#			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperNu,colour=as.factor(fillColor)),size=1.5)
#			plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
#				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
#				theme(legend.position="none")+
#				labs(x="Cluster",y="Shape Parameter")
#			plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
#			plotObj<-plotObj+labs(title="",plot.title=element_text(size=10))
#			# Margin order is (top,right,bottom,left)
#			plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
#			print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=2))
		} else {
			rownames(riskDF)<-seq(1,nrow(riskDF),by=1)

			# Create the risk plot
			plotObj<-ggplot(riskDF)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=risk,yintercept=meanRisk))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=risk,fill=as.factor(fillColor)),outlier.size=0.5)
			if (!is.null(riskLim)) plotObj<-plotObj+coord_cartesian(ylim = riskLim)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerRisk,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperRisk,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+
				labs(x="Cluster",y=ifelse(showRelativeRisk,'RR',
				ifelse(yModel=="Categorical"||yModel=="Bernoulli"||yModel=="Binomial","Probability","E[Y]")))
			plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
			plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
			# Margin order is (top,right,bottom,left)
			plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=2))
		}	
	}

	# Create a bar chart of cluster empiricals
	if((!is.null(yModel))){
		if(yModel!="Categorical"){
			plotObj<-ggplot(empiricalDF)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=empiricals,colour=as.factor(fillColor)),size=3)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=empiricals,yintercept=meanEmpirical))
			plotObj<-plotObj+scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")
			plotObj<-plotObj+labs(title='Empirical Data',plot.title=element_text(size=10))
			plotObj<-plotObj+theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10,angle=90))
			plotObj<-plotObj+
				labs(y=ifelse(yModel=="Bernoulli","Proportion of cases",
				ifelse(yModel=="Binomial","Avg Proportion of occurrence",
				ifelse(yModel=="Poisson","Avg Count",
				ifelse(yModel=="Survival","Avg Survival Time",
				ifelse(yModel=="Categorical","Avg Proportion of occurrence","Avg Y"))))),x="Cluster")
			plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.15,0.5,0.5,1),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=1))
		}
	}
	# Create a bar chart of cluster sizes
	plotObj<-ggplot(sizeDF)
	plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=clusterSize,colour=as.factor(fillColor)),size=3)
	plotObj<-plotObj+scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+theme(legend.position="none")
	plotObj<-plotObj+labs(title="Size",plot.title=element_text(size=10))
	plotObj<-plotObj+theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10,angle=90))
	plotObj<-plotObj+labs(y="No. of Subjects",x="Cluster")
	plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.15,0.5,0.5,1),'lines'))
	print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=1))

	# Loop over the covariates
	for(j in 1:nCovariates){
		if(xModel=='Discrete'){
			profileDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
				"lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
			for(k in 1:nCategories[j]){
				probMat<-profile[,meanSortIndex,j,k]
				nPoints<-nrow(probMat)
				probMeans<-apply(probMat,2,mean)
				probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
				probLower<-apply(probMat,2,quantile,0.05)
				probUpper<-apply(probMat,2,quantile,0.95)
		
				# Get the plot colors
				probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
				ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))
			
	
				for(c in whichClusters){
					profileDF<-rbind(profileDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
						"category"=rep(k-1,nPoints),
						"meanProb"=rep(probMean,nPoints),
						"lowerProb"=rep(probLower[c],nPoints),
						"upperProb"=rep(probUpper[c],nPoints),
						"fillColor"=rep(probColor[c],nPoints)))
						 rownames(profileDF)<-seq(1,nrow(profileDF),1)
				
				}
			}
		
			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
			plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
			plotObj<-plotObj+theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
				theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=j+2))
		}else if(xModel=='Normal'){
			# Plot the means
			profileDF<-data.frame("mu"=c(),"cluster"=c(),"muMean"=c(),
				"lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
			muMat<-profile[,meanSortIndex,j]
			muMeans<-apply(muMat,2,mean)
			muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
			muLower<-apply(muMat,2,quantile,0.05)
			muUpper<-apply(muMat,2,quantile,0.95)
			# The next line is to avoid outliers spoiling plot scales
			plotMax<-max(muUpper)
			plotMin<-min(muLower)
	
			# Get the plot colors
			muColor<-ifelse(muLower>rep(muMean,length(muLower)),"high",
			ifelse(muUpper<rep(muMean,length(muUpper)),"low","avg"))
			for(c in whichClusters){
				plotMu<-muMat[,c]
				plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
				nPoints<-length(plotMu)
				profileDF<-rbind(profileDF,data.frame("mu"=plotMu,"cluster"=rep(c,nPoints),
					"meanMu"=rep(muMean,nPoints),
					"lowerMu"=rep(muLower[c],nPoints),
					"upperMu"=rep(muUpper[c],nPoints),
					"fillColor"=rep(muColor[c],nPoints)))
			}

			rownames(profileDF)<-seq(1,nrow(profileDF),1)
			
			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
				plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
				plotObj<-plotObj+
					theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
					theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j+2))
		
			# Plot the variances
			profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
				"lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
			sigmaMat<-profileStdDev[,meanSortIndex,j,j]
			sigmaMeans<-apply(sigmaMat,2,mean)
			sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
			sigmaLower<-apply(sigmaMat,2,quantile,0.05)
			sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
			# The next line is to avoid outliers spoiling plot scales
			plotMax<-max(sigmaUpper)
	
			# Get the plot colors
			sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,length(sigmaLower)),"high",
			ifelse(sigmaUpper<rep(sigmaMean,length(sigmaUpper)),"low","avg"))
			for(c in whichClusters){
				plotSigma<-sigmaMat[,c]
				plotSigma<-plotSigma[plotSigma<plotMax]
				nPoints<-length(plotSigma)
				profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
					"meanSigma"=rep(sigmaMean,nPoints),
					"lowerSigma"=rep(sigmaLower[c],nPoints),
					"upperSigma"=rep(sigmaUpper[c],nPoints),
					"fillColor"=rep(sigmaColor[c],nPoints)))
			}
			rownames(profileDF)<-seq(1,nrow(profileDF),1)

			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
			plotObj<-plotObj+
				theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
				theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j+2))
		}else if(xModel=='Mixed'){
			if (j<=nDiscreteCovs){
			profileDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
				"lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
			for(k in 1:nCategories[j]){
				if (nDiscreteCovs==1) {
					probMat<-profilePhi[,meanSortIndex,1,k]
				} else {
					probMat<-profilePhi[,meanSortIndex,j,k]
				}
				nPoints<-nrow(probMat)
				probMeans<-apply(probMat,2,mean)
				probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
				probLower<-apply(probMat,2,quantile,0.05)
				probUpper<-apply(probMat,2,quantile,0.95)
		
				# Get the plot colors
				probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
				ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))
			
	
				for(c in whichClusters){
					profileDF<-rbind(profileDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
						"category"=rep(k-1,nPoints),
						"meanProb"=rep(probMean,nPoints),
						"lowerProb"=rep(probLower[c],nPoints),
						"upperProb"=rep(probUpper[c],nPoints),
						"fillColor"=rep(probColor[c],nPoints)))
						 rownames(profileDF)<-seq(1,nrow(profileDF),1)
				
				}
			}
		
			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
			plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
			plotObj<-plotObj+theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
				theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=j+2))
			} else {			
			# Plot the means
			profileDF<-data.frame("mu"=c(),"cluster"=c(),"muMean"=c(),
				"lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
			if (nContinuousCovs==1){
				muMat<-profileMu[,meanSortIndex,1]
			} else {
				muMat<-profileMu[,meanSortIndex,(j-nDiscreteCovs)]
			}
			muMeans<-apply(muMat,2,mean)
			muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
			muLower<-apply(muMat,2,quantile,0.05)
			muUpper<-apply(muMat,2,quantile,0.95)
			# The next line is to avoid outliers spoiling plot scales
			plotMax<-max(muUpper)
			plotMin<-min(muLower)
			
			# Get the plot colors
			muColor<-ifelse(muLower>rep(muMean,length(muLower)),"high",
			ifelse(muUpper<rep(muMean,length(muUpper)),"low","avg"))
			for(c in whichClusters){
				plotMu<-muMat[,c]
				plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
				nPoints<-length(plotMu)
				profileDF<-rbind(profileDF,data.frame("mu"=plotMu,"cluster"=rep(c,nPoints),
					"meanMu"=rep(muMean,nPoints),
					"lowerMu"=rep(muLower[c],nPoints),
					"upperMu"=rep(muUpper[c],nPoints),
					"fillColor"=rep(muColor[c],nPoints)))
			}
			rownames(profileDF)<-seq(1,nrow(profileDF),1)
			
			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
				plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
				plotObj<-plotObj+
					theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
					theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j+2))
		
			# Plot the variances
			profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
				"lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
			if (nContinuousCovs==1){
				sigmaMat<-profileStdDev[,meanSortIndex,1,1]
			} else {
				sigmaMat<-profileStdDev[,meanSortIndex,(j-nDiscreteCovs),(j-nDiscreteCovs)]
			}
			sigmaMeans<-apply(sigmaMat,2,mean)
			sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
			sigmaLower<-apply(sigmaMat,2,quantile,0.05)
			sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
			# The next line is to avoid outliers spoiling plot scales
			plotMax<-max(sigmaUpper)
	
			# Get the plot colors
			sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,length(sigmaLower)),"high",
			ifelse(sigmaUpper<rep(sigmaMean,length(sigmaUpper)),"low","avg"))
			for(c in whichClusters){
				plotSigma<-sigmaMat[,c]
				plotSigma<-plotSigma[plotSigma<plotMax]
				nPoints<-length(plotSigma)
				profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
					"meanSigma"=rep(sigmaMean,nPoints),
					"lowerSigma"=rep(sigmaLower[c],nPoints),
					"upperSigma"=rep(sigmaUpper[c],nPoints),
					"fillColor"=rep(sigmaColor[c],nPoints)))
			}
			rownames(profileDF)<-seq(1,nrow(profileDF),1)
			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
			plotObj<-plotObj+
				theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
				theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j+2))
		}
	
		}
	}
	dev.off()
	return(meanSortIndex)
}
	

# Calculate predictions, and if possible assess predictive performance
calcPredictions<-function(riskProfObj,predictResponseFileName=NULL, doRaoBlackwell=F, fullSweepPredictions=F,fullSweepLogOR=F,fullSweepHazardRatio=F){

	riskProfClusObj=NULL
	clusObjRunInfoObj=NULL
	yModel=NULL
	reportBurnIn=NULL
	nBurn=NULL
	nFilter=NULL
	nSweeps=NULL
	nPredictSubjects=NULL
	fullPredictFile=NULL
	nFixedEffects=NULL
	directoryPath=NULL
	fileStem=NULL
	nCategoriesY=NULL
	nSubjects=NULL
	weibullFixedShape=NULL
	
	for (i in 1:length(riskProfObj)) assign(names(riskProfObj)[i],riskProfObj[[i]])
	for (i in 1:length(riskProfClusObj)) assign(names(riskProfClusObj)[i],riskProfClusObj[[i]])
	for (i in 1:length(clusObjRunInfoObj)) assign(names(clusObjRunInfoObj)[i],clusObjRunInfoObj[[i]])
	
	if(yModel=="Poisson"||yModel=="Normal"){
		if (fullSweepLogOR==T){
			fullSweepLogOR=F
			cat("Log odds ratio does not make sense for Poisson or Normal response\n")
		}
	}


	if(!yModel=="Survival"){
		if (fullSweepHazardRatio==T){
			fullSweepHazardRatio=F
			cat("Hazard ratio makes sense for Survival response only\n")
		}
		if (!weibullFixedShape){
			doRaoBlackwell=F
			cat("For Survival response with cluster specific shape parameter Rao-Blackwell predictions are not possible.\n")
		}
	}
	
	if (yModel=="Survival"){
		restrictedMeanSurvival<-max(clusObjRunInfoObj$yMat[,1])
	} 

	firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
	lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter
	nSamples<-lastLine-firstLine+1
	
	# First of all we see if there a covariate file has been supplied
	if(nPredictSubjects==0){
		stop("No prediction subjects processed by C++\n")
	}
	
	# Get the response and fixed effects data if available
	responseProvided<-F
	extraInfoProvided<-F # extra info is denominator for binomial and poisson
	fixedEffectsProvided<-F
	if (is.null(predictResponseFileName)&&fullPredictFile==TRUE){
		predictResponseFileName = file.path(directoryPath,paste(fileStem,'_predictFull.txt',sep=''))
	}
	if(!is.null(predictResponseFileName)){
		predictResponseData<-scan(predictResponseFileName,quiet=T)
		predictResponseMat<-matrix(predictResponseData[2:length(predictResponseData)],
			nrow=nPredictSubjects,byrow=T)
		predictYMat<-matrix(predictResponseMat[,1],ncol=1)
		if(yModel=="Poisson"||yModel=="Binomial"||yModel=="Survival"){
			predictYMat<-cbind(predictYMat,predictResponseMat[,ncol(predictResponseMat)])
			if(all(predictYMat[,2]>-999)){
				extraInfoProvided<-T
			}
		}
		if(all(predictYMat[,1]>-999)){
			responseProvided<-T
		}
		if(nFixedEffects>0){
			fixedEffectsProvided<-T
			predictWMat<-matrix(predictResponseMat[,2:(nFixedEffects+1)],nrow=nPredictSubjects)
			# Set missing values to their average or reference value
			predictWMat[predictWMat==-999]<-0
		}
	}

	if(fixedEffectsProvided){
		betaFileName <-file.path(directoryPath,paste(fileStem,'_beta.txt',sep=''))
		betaFile<-file(betaFileName,open="r")
		betaArray<-array(0,dim=c(nSamples,nFixedEffects,nCategoriesY))
		for(sweep in firstLine:lastLine){
			if(sweep==firstLine){
				skipVal<-firstLine-1
			}else{
				skipVal<-0
			}
			if (yModel=="Categorical") {
				currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects*(nCategoriesY-1),quiet=T)
				currBeta<-matrix(currBetaVector,ncol=(nCategoriesY-1),byrow=T)
				currBeta<-cbind(rep(0,dim(currBeta)[1]),currBeta)
			} else {
				currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects,quiet=T)
				currBeta<-matrix(currBetaVector,ncol=nCategoriesY,byrow=T)
			}
			betaArray[sweep-firstLine+1,,]<-currBeta
		}
	}

	if (yModel=="Survival"){
		if (weibullFixedShape){
			nuFileName<-file.path(directoryPath,paste(fileStem,'_nu.txt',sep=''))
			nuFile<-file(nuFileName,open="r")
			nu<-mean(read.table(nuFile)[,1])
			close(nuFile)
		}
	} 

	# Already done the allocation in the C++
	if(doRaoBlackwell){
		# Construct the RB theta file name
		thetaFileName<-file.path(directoryPath,paste(fileStem,'_predictThetaRaoBlackwell.txt',sep=''))
		thetaArray<-array(0,dim=c(nSamples,nPredictSubjects,nCategoriesY))
		# Read the RB theta data
		if (yModel=="Categorical"){
			thetaMat<-matrix(scan(thetaFileName,what=double(),quiet=T),byrow=T,ncol=nPredictSubjects*(nCategoriesY-1))
			for (sweep in firstLine:lastLine){
				thetaArrayRow<-cbind(rep(0,nPredictSubjects),matrix(thetaMat[sweep,],ncol=nCategoriesY-1,byrow=T))
				thetaArray[sweep-firstLine+1,,]<-thetaArrayRow
			}
		} else {
			thetaMat<-matrix(scan(thetaFileName,what=double(),quiet=T),byrow=T,ncol=nPredictSubjects)
			thetaArray[,,1]<-thetaMat[firstLine:nrow(thetaMat),]
		}
	}else{
		# Construct the file names
		zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
		zFile<-file(zFileName,open="r")
		thetaFileName<-file.path(directoryPath,paste(fileStem,'_theta.txt',sep=''))
		thetaFile<-file(thetaFileName,open="r")
		nClustersFileName<-file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
		nClustersFile<-file(nClustersFileName,open="r")
		if (yModel=="Survival"&&!weibullFixedShape){
			nuFileName<-file.path(directoryPath,paste(fileStem,'_nu.txt',sep=''))
			nuFile<-file(nuFileName,open="r")
			nuArrayPred<-array(0,dim=c(nSamples,nPredictSubjects))
		}
		# initialise theta and beta arrays
		thetaArray<-array(0,dim=c(nSamples,nPredictSubjects,nCategoriesY))
		for(sweep in firstLine:lastLine){
			if(sweep==firstLine){
				skipVal<-firstLine-1
			}else{
				skipVal<-0
			}
			currZ<-1+scan(zFile,what=integer(),skip=skipVal,n=nSubjects+nPredictSubjects,quiet=T)
			tmpCurrZ<-currZ
			currZ<-currZ[(nSubjects+1):length(currZ)]
			currNClusters<-scan(nClustersFile,skip=skipVal,n=1,what=integer(),quiet=T)
			if (yModel=="Categorical") {
				currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currNClusters*(nCategoriesY-1),quiet=T)
				currTheta<-matrix(currThetaVector,ncol=(nCategoriesY-1),byrow=T)
				currTheta<-cbind(rep(0,dim(currTheta)[1]),currTheta)
			} else {
				currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currNClusters,quiet=T)
				currTheta<-matrix(currThetaVector,ncol=nCategoriesY,byrow=T)
				if (yModel=="Survival"&&!weibullFixedShape) currNu<-scan(nuFile,what=double(),skip=skipVal,n=currNClusters,quiet=T)
			}
			thetaArray[sweep-firstLine+1,,]<-as.matrix(currTheta[currZ,],ncol=nCategoriesY)
			if (yModel=="Survival"&&!weibullFixedShape) nuArrayPred[sweep-firstLine+1,]<-as.vector(currNu[currZ])

		}
		close(zFile)
		close(thetaFile)
		close(nClustersFile)
		if (yModel=="Survival"&&!weibullFixedShape) close(nuFile)
	}


	
	# Use the values of theta to derive the predicted values
	predictedY<-array(0,dim=c(nSamples,nPredictSubjects,nCategoriesY))
	for(sweep in 1:nSamples){
		if(sweep==1||sweep%%1000==0){
			cat(paste("Processing sweep",sweep,"of ",nSamples,"\n"))
		}
	
		lambda<-thetaArray[sweep,,]
		if(fixedEffectsProvided){
			currBeta<-betaArray[sweep,,]
			lambda<-lambda+predictWMat%*%currBeta
		}
	
		if(yModel=='Poisson'){
			if(extraInfoProvided){
				# Add in the offset
				lambda<-lambda+log(predictYMat[,2])
			}
		}
	
		if(yModel=='Bernoulli'){
			predictedY[sweep,,]<-exp(lambda)/(1+exp(lambda))
		}else if(yModel=='Binomial'){
			if(extraInfoProvided){
				predictedY[sweep,,]<-predictYMat[,2]*exp(lambda)/(1+exp(lambda))
			}else{
				predictedY[sweep,,]<-exp(lambda)/(1+exp(lambda))
			}
		}else if(yModel=='Poisson'){
			predictedY[sweep,,]<-exp(lambda)
		}else if(yModel=='Normal'){
			predictedY[sweep,,]<-lambda
		}else if(yModel=='Survival'){
			if (!weibullFixedShape) nu<-nuArrayPred[sweep,]
			tmpPredictedY<-1/((exp(lambda))^(1/nu))*gamma(1+1/nu)
			tmpPredictedY[which(tmpPredictedY>restrictedMeanSurvival)]<-restrictedMeanSurvival
			predictedY[sweep,,]<-tmpPredictedY
		}else if(yModel=="Categorical"){
			predictedY[sweep,,]<-exp(lambda)/rowSums(exp(lambda))
		}
	}

	if(responseProvided){
		bias<-apply(predictedY,2,median)-predictYMat[,1]
		rmse<-sqrt(mean(bias^2))
		mae<-mean(abs(bias))
		bias<-mean(bias)
		output<-list("bias"=bias,"rmse"=rmse,
			"observedY"=predictYMat[,1],
			"predictedY"=apply(predictedY,c(2,3),median),
			"doRaoBlackwell"=doRaoBlackwell,
			"mae"=mae)
	}else{
		output<-list("bias"=NA,"rmse"=NA,"observedY"=NA,"predictedY"=apply(predictedY,c(2,3),median),"doRaoBlackwell"=doRaoBlackwell,
			"mae"=NA)
	}
	if(fullSweepPredictions){
		output$predictedYPerSweep<-predictedY
	}
	if(fullSweepLogOR){
		logORPerSweep<-matrix(0,ncol=ncol(predictedY),nrow=nrow(predictedY))
		for(i in 1:dim(thetaArray)[1]){
			# Always relative to the first prediction subject
			logORPerSweep[i,]<-thetaArray[i,,]-thetaArray[i,1,]
		}
		output$logORPerSweep<-logORPerSweep
	}
	if (fullSweepHazardRatio){
		fullHR<-matrix(0,ncol=ncol(predictedY),nrow=nrow(predictedY))
		for(i in 1:dim(thetaArray)[1]){
			# Always relative to the first prediction subject
			fullHR[i,]<-predictedY[i,,]/predictedY[i,1,]
		}
		output$fullHR<-fullHR

	}
	if(fixedEffectsProvided){
		close(betaFile)
	}
	return(output)
}

# Show the continuous hyperparameter for variable selection
summariseVarSelectRho<-function(runInfoObj){
	
	directoryPath=NULL
	fileStem=NULL
	nCovariates=NULL
	reportBurnIn=NULL
	nBurn=NULL
	nFilter=NULL
	nSweeps=NULL

	for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])

	# Rho file name
	rhoFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
	rhoMat<-matrix(scan(rhoFileName,what=double(),quiet=T),ncol=nCovariates,byrow=T)
	
	# Restrict to after burn in
	firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
	lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter

	rhoMat<-rhoMat[firstLine:lastLine,]
	
	rhoMean<-apply(rhoMat,2,mean)
	rhoMedian<-apply(rhoMat,2,median)
	rhoLowerCI<-apply(rhoMat,2,quantile,0.05)
	rhoUpperCI<-apply(rhoMat,2,quantile,0.95)

	return(list("rho"=rhoMat,"rhoMean"=rhoMean,"rhoMedian"=rhoMedian,"rhoLowerCI"=rhoLowerCI,"rhoUpperCI"=rhoUpperCI))
	
}
	

# Function to compute the marginal model posterior (only for discrete covariates and Bernoulli outcome)
margModelPosterior<-function(runInfoObj,allocation){

	xModel=NULL
	yModel=NULL
	varSelect=NULL
	nSubjects=NULL
	xMat=NULL
	directoryPath=NULL
	fileStem=NULL
	includeResponse=NULL
	nFixedEffects=NULL
	reportBurnIn=NULL
	nBurn=NULL
	nFilter=NULL
	nSweeps=NULL
	nProgress=NULL
	dPitmanYor=NULL
	nCovariates=NULL

	for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])

	# this function only works for Bernoulli outcome and discrete covariates, so check that it is used correctly
	if (includeResponse) {
		if (xModel!="Discrete"||yModel!="Bernoulli") stop("ERROR: The computation of the marginal model posterior has only been implemented for Bernoulli outcome and discrete covariates.")
	} else {
		if (xModel!="Discrete") stop("ERROR: The computation of the marginal model posterior has only been implemented for Bernoulli outcome and discrete covariates.")
	}
	# no variable selection has been implemented
	if (varSelect==TRUE) print("Warning: Variable selection is not taken into account for the computation of marginal model posterior")
	# the subjects with missing values are simply removed for now
	missingX<-FALSE
	for (k in 1:nSubjects){
		if (nCovariates>1&&sum(xMat[k,]==-999)>0) {
			missingX<-TRUE
			stop("ERROR: No missing value handling technique has been implemented for the marginal model posterior. This function cannot be run if missing values are present.")
		}
		if (nCovariates==1&&sum(xMat[k]==-999)>0) {
			missingX<-TRUE
			stop("ERROR: No missing value handling technique has been implemented for the marginal model posterior. This function cannot be run if missing values are present.")
		}
	}
	# this function only works for dPitmanYor == 0, ie only for Dirichlet process prior
	if (runInfoObj$dPitmanYor !=0) stop("ERROR: the marginal model posterior has only been implemented for Dirichlet process priors.")

	# read in value of hyperparameters
	runData<-readLines(file.path(directoryPath,paste(fileStem,'_log.txt',sep='')))

	hyperParams<-list()

	if (includeResponse==T){
		sigmaTheta<-runData[grep('sigmaTheta',runData)]
		sigmaTheta<-substr(sigmaTheta,regexpr(':',sigmaTheta)+1,nchar(sigmaTheta))
		sigmaTheta<-gsub(' ','',sigmaTheta)
		sigmaTheta<-gsub('\t','',sigmaTheta)
		sigmaTheta<-as.integer(sigmaTheta)

		dofTheta<-runData[grep('dofTheta',runData)]
		dofTheta<-substr(dofTheta,regexpr(':',dofTheta)+1,nchar(dofTheta))
		dofTheta<-gsub(' ','',dofTheta)
		dofTheta<-gsub('\t','',dofTheta)
		dofTheta<-as.integer(dofTheta)

		hyperParams<-list(sigmaTheta=sigmaTheta,dofTheta=dofTheta)
	}
	
	if (nFixedEffects>0){
		sigmaBeta<-runData[grep('sigmaBeta',runData)]
		sigmaBeta<-substr(sigmaBeta,regexpr(':',sigmaBeta)+1,nchar(sigmaBeta))
		sigmaBeta<-gsub(' ','',sigmaBeta)
		sigmaBeta<-gsub('\t','',sigmaBeta)
		sigmaBeta<-as.integer(sigmaBeta)

		dofBeta<-runData[grep('dofBeta',runData)]
		dofBeta<-substr(dofBeta,regexpr(':',dofBeta)+1,nchar(dofBeta))
		dofBeta<-gsub(' ','',dofBeta)
		dofBeta<-gsub('\t','',dofBeta)
		dofBeta<-as.integer(dofBeta)

		hyperParams$sigmaBeta<-sigmaBeta
		hyperParams$dofBeta<-dofBeta
	}


	# set alpha, whether it has been estimated or it is fixed
	# find out if it was fixed or estimated 
	alphaUpdate<-runData[grep('Update alpha',runData)]
	alphaUpdate<-substr(alphaUpdate,regexpr(':',alphaUpdate)+1,nchar(alphaUpdate))
	alphaUpdate<-gsub(' ','',alphaUpdate)
	if (alphaUpdate=="False") {
		alpha<-runData[grep('Fixed alpha: ',runData)]
		alpha<-substr(alpha,regexpr(':',alpha)+1,nchar(alpha))
		alpha<-as.integer(alpha)
	} else {
		# if alpha wasn't fixed, take median value of chain
		firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
		skipLines<-ifelse(reportBurnIn,nBurn/nFilter+1,0)
		lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter		
		alphaFileName <- file(file.path(directoryPath,paste(fileStem,'_alpha.txt',sep='')))
		open(alphaFileName)
		alphaValues<-vector()
		alphaValues[1]<-scan(alphaFileName,what=double(),skip=skipLines,nlines=1,quiet=T)
		for (i in (firstLine+1):lastLine){
			alphaValues[i-firstLine]<-scan(alphaFileName,what=double(),skip=0,nlines=1,quiet=T)
		}
		close(alphaFileName)
		alpha<-median(alphaValues)
	}
	runInfoObj$alphaMPP <- alpha

	if (xModel=="Discrete"){
		aPhi<-runData[grep('aPhi',runData)]
		aPhi<-substr(aPhi,regexpr(':',aPhi)+1,nchar(aPhi))
		aPhi<-strsplit(aPhi," ")[[1]][-1]
		aPhi<-as.integer(aPhi)
		hyperParams$aPhi<-aPhi
	}

	runInfoObj$hyperParams <- hyperParams

	# read first allocation iteration after burnin
	firstLine<-ifelse(reportBurnIn,nBurn/nFilter+1,1)
	skipLines<-ifelse(reportBurnIn,nBurn/nFilter,0)
	lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter	

	# open allocation file
	if (missing(allocation)){
		zFileName <- file(file.path(directoryPath,paste(fileStem,'_z.txt',sep='')))
		open(zFileName)
		zAllocCurrent<-scan(zFileName,what=integer(),skip=skipLines,nlines=1,quiet=T)
		zAllocCurrent<-zAllocCurrent[1:nSubjects]
	} else {
		zAllocCurrent<-allocation[1:nSubjects]
		firstLine<-0
		skipLines<-0
		lastLine<-0
	}
	

	# initialise output vectors
	margModPost<-rep(0,length=(lastLine-firstLine+1))

	clusterSizes<-table(zAllocCurrent)
	nClusters<-length(clusterSizes)
	# set initial values of parameters of interest
	if (nFixedEffects>0){
		parFirstIter<-c(rep(0,nClusters),rep(0,nFixedEffects))
	} else {
		parFirstIter<-c(rep(0,nClusters))
	}
	# compute marginal model posterior

	output<-.pZpXpY(zAlloc=zAllocCurrent, par=parFirstIter, clusterSizes=clusterSizes, nClusters=nClusters, runInfoObj=runInfoObj, alpha=alpha)
	margModPost[1]<-output$margModPost
	if (missing(allocation)){
		for (iter in (firstLine+1):lastLine){
			if (iter%%nProgress==0) print(iter)
			# identify allocations for this sweep
			zAllocCurrent<-scan(zFileName,what=integer(),nlines=1,quiet=T)
			zAllocCurrent<-zAllocCurrent[1:nSubjects]
			# parameters
			# number of elements in each cluster
			clusterSizes<-table(zAllocCurrent)
			# number of clusters
			nClusters<-length(clusterSizes)
			# computing the marginal likelihood
			# version using the previous beta mode for next step	
			if (nFixedEffects>0){
				parTmp<-c(rep(0,nClusters),rep(0,nFixedEffects))
			} else {
				parTmp<-c(rep(0,nClusters))
			}
	
			output<-.pZpXpY(zAlloc=zAllocCurrent,par=parTmp, clusterSizes=clusterSizes, nClusters=nClusters, runInfoObj = runInfoObj, alpha=alpha)
			margModPost[iter-firstLine+1]<-output$margModPost

		}	
	}
	if (missing(allocation)){
		close(zFileName)
	}

	write.table(margModPost,file.path(directoryPath,paste(fileStem,"_margModPost.txt",sep="")), col.names = FALSE,row.names = FALSE)
	return(list("meanMargModPost"=mean(margModPost),"runInfoObj"=runInfoObj))
}

# internal function
# Function to evaluate pYGivenZW for Bernoulli outcome - required for the Laplace approximation
.pYGivenZW_Bernoulli<-function(thetaBeta,zAlloc,nClusters,hyperParams,
		yMat,wMat,nSubjects,nFixedEffects,nTableNames,constants){
	dimThetaBeta<-length(thetaBeta)
	theta<-head(thetaBeta,nClusters)
	if (nFixedEffects>0){
		beta<-tail(thetaBeta,-nClusters)
		betaW<-as.matrix(wMat)%*%beta
	} else {
		beta<-0
		hyperParams$dofBeta<-0
		hyperParams$sigmaBeta<-0
		betaW<-0
	}
	nTableNames<-as.integer(nTableNames)
	maxNTableNames<-max(nTableNames)
	out<-.Call('pYGivenZW',beta,theta,zAlloc,hyperParams$sigmaBeta,
		hyperParams$sigmaTheta,hyperParams$dofTheta,hyperParams$dofBeta,nSubjects,
		yMat,betaW,nFixedEffects,nTableNames,constants,
		maxNTableNames,PACKAGE = 'PReMiuM')
	return(out)
}

# internal function
# Function to evaluate the gradient of pYGivenZW for Bernoulli outcome - required for the Laplace approximation
.pYGivenZW_Bernoulli_gradient<-function(thetaBeta,zAlloc,nClusters,hyperParams,
		yMat,wMat,nSubjects,nFixedEffects,nTableNames,constants){
	dimThetaBeta<-length(thetaBeta)
	theta<-head(thetaBeta,nClusters)
	nTableNames<-as.integer(nTableNames)
	maxNTableNames<-max(nTableNames)
	if (nFixedEffects>0){
		beta<-tail(thetaBeta,-nClusters)
		betaW<-as.matrix(wMat)%*%beta
		yPred<-.Call('GradpYGivenZW',beta,theta,zAlloc,nSubjects,
			betaW,yMat,nFixedEffects,nTableNames,
			maxNTableNames,PACKAGE = 'PReMiuM')
		wPred<- as.matrix(wMat)*yPred
	} else {
		beta<-0
		betaW<-0
		yPred<-.Call('GradpYGivenZW',beta,theta,zAlloc,nSubjects,
			betaW,yMat,nFixedEffects,nTableNames,
			maxNTableNames,PACKAGE = 'PReMiuM')
	}
	gr.theta<-rep(0,nClusters)
	for (k in 1:nClusters){
		gr.theta[k]<-sum(yPred[zAlloc==nTableNames[k]])
	}
	# contribution from the derivative of theta
	numeratorTheta<-(hyperParams$dofTheta+1)*theta
	paramsTheta<-hyperParams$dofTheta*hyperParams$sigmaTheta^2
	gr.theta<-gr.theta-numeratorTheta/(paramsTheta+theta^2)
	out<- -c(gr.theta)/nSubjects
	if (nFixedEffects>0){
	# contribution from the derivative of beta		
		gr.beta<-rep(0,nFixedEffects)
		gr.beta<-gr.beta+apply(wPred,2,sum)
		numeratorBeta<-(hyperParams$dofBeta+1)*beta
		paramsBeta<-hyperParams$dofBeta*hyperParams$sigmaBeta^2
		gr.beta<-gr.beta-numeratorBeta/(paramsBeta+beta^2)
		out<- c(out,-gr.beta/nSubjects)	
	}
	return(out)
}

# internal function
# Function to evaluate the Hessian matrix of pYGivenZW for Bernoulli outcome - required for the Laplace approximation
.pYGivenZW_Bernoulli_hessianMat<-function(thetaBeta,zAlloc,nClusters,hyperParams,
		wMat,nSubjects,nFixedEffects,nSweeps,nTableNames){
	dimThetaBeta<-length(thetaBeta)
	theta<-head(thetaBeta,nClusters)
	if (nFixedEffects>0){
		wMat<-as.matrix(wMat)
		beta<-tail(thetaBeta,-nClusters)
		betaW<-wMat%*%beta
	} else {
		beta<-0
		betaW<-0
	}
	hes.mat<-matrix(0,ncol=dimThetaBeta,nrow=dimThetaBeta)
	exp.predictor<-rep(0,nSubjects)
	denomPred<-rep(0,nSubjects)
	predictor<-rep(0,nSubjects)
	thetaTmp<-rep(0,length=max(as.integer(nTableNames)))
	k<-1
	for (i in as.integer(nTableNames)+1) {
		thetaTmp[i]<-theta[k]	
		k<-k+1
	}
	for (i in 1:nSubjects){
		predictor[i]<-thetaTmp[zAlloc[i]+1]
	}
	if (nFixedEffects>0){
		exp.predictor<-exp(predictor+as.vector(betaW))
	} else {
		exp.predictor<-exp(predictor)
	}
	denomPred<-exp.predictor/(exp.predictor+1)^2
	# second derivatives
	# wrt theta_i, theta_j is = 0
	# wrt theta_i, theta_i is
	paramsTheta<-hyperParams$dofTheta*hyperParams$sigmaTheta^2
	for (k in 1:nClusters){
		diag(hes.mat)[k]<-sum(denomPred[zAlloc==nTableNames[k]]) + (hyperParams$dofTheta+1)*
			(paramsTheta - theta[k]^2)/
			(paramsTheta + theta[k]^2)^2
     
	}
	if (nFixedEffects>0){
		# wrt beta_i, theta_j is
		for (k in 1:nFixedEffects){
			for (j in 1:nClusters){
				indexJ<-zAlloc==nTableNames[j]
				hes.mat[j,(nClusters+k)]<-hes.mat[(nClusters+k),j]<-
					sum(wMat[indexJ,k]*denomPred[indexJ])
			}
		}
		# wrt beta_i, beta_j is
		if (nFixedEffects>1){
			for (k in 1:(nFixedEffects-1)){
				for (j in (k+1):nFixedEffects){
					hes.mat[(nClusters+k),(nClusters+j)]<-hes.mat[(nClusters+j),(nClusters+k)]<-
						sum((wMat[,k]*wMat[,j])*denomPred)
				}
			}
		}
		# wrt beta_i, beta_i is
		betaSq<-beta*beta
		paramsBeta<-hyperParams$dofBeta*hyperParams$sigmaBeta^2
		diag(hes.mat)[(nClusters+1):(dimThetaBeta)]<- 
			apply(wMat^2 * denomPred,2,sum)+(hyperParams$dofBeta+1)*
			(paramsBeta - betaSq)/
			(paramsBeta + betaSq)^2
	}
	return(hes.mat/nSubjects) 
}	


# internal function
# marginal model posterior for one iteration
.pZpXpY<-function(zAlloc,runInfoObj,par=NULL,clusterSizes,nClusters,alpha){

	nSweeps=NULL
	nFixedEffects=NULL
	nCategories=NULL
	hyperParams=NULL
	nCovariates=NULL
	xMat=NULL
	includeResponse=NULL
	yMat=NULL
	wMat=NULL
	nSubjects=NULL
	alphaMPP=NULL

	for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])

	# marginal likelihood is p(Z|alpha,X,Y)
	
	# p(Z|alpha,X,Y) is prop to p(X|Z)p(Y|Z,W)p(Z|alpha)
	# we are going to call the three quantities pX, pY and pZ
	
	nTableNames<-names(clusterSizes)
	
	# set initial parameters at 0 if a better starting point is not provided	
	if (is.null(par)) par<-c(rep(0,nClusters),rep(0,nFixedEffects))

	# number of elements in following clusters
	nPlus<-head(c(rev(cumsum(rev(clusterSizes[-1]))),0),-1)	

	# computation of pX+pZ
	pZpX<-.Call('pZpX',nClusters,nCategories,hyperParams$aPhi,clusterSizes,nCovariates, zAlloc, as.vector(as.matrix(xMat)), as.integer(nTableNames), alphaMPP, nPlus, PACKAGE = 'PReMiuM')

	# computation of pY
	if (includeResponse==T){
		# constants of pY
		constants<-nClusters*(0.5*(hyperParams$dofTheta+1)*log(hyperParams$dofTheta)-
			lgamma(0.5*(hyperParams$dofTheta+1))+0.5*log(hyperParams$dofTheta*pi)+lgamma(hyperParams$dofTheta*0.5))
		if (nFixedEffects>0) constants<-constants+nFixedEffects*(0.5*(hyperParams$dofBeta+1)*log(hyperParams$dofBeta)-
			lgamma(0.5*(hyperParams$dofBeta+1))+0.5*log(hyperParams$dofBeta*pi)+lgamma(hyperParams$dofBeta*0.5))
		# computation of min for Laplace approximation		
		laplaceOut <- optim(par = par, # initial parameters
			fn = .pYGivenZW_Bernoulli, # function to minimise over
			control=list(fnscale=1), # need to minimise 
			method="L-BFGS-B", # fastest method according to tests
			gr=.pYGivenZW_Bernoulli_gradient, # gradient
			hessian=FALSE, # we don't require the approximate computation of the hessian matrix, it's computed separately and exactly
			# other parameters necessary for fn above
			nClusters=nClusters, 
			zAlloc=zAlloc,hyperParams=hyperParams,
			yMat=yMat,wMat=wMat,nSubjects=nSubjects,
			nFixedEffects=nFixedEffects,nTableNames=nTableNames,constants=constants)
		laplaceMode <- laplaceOut$value 
		# computation of the hessian matrix		
		laplaceHessian <- .pYGivenZW_Bernoulli_hessianMat(laplaceOut$par,zAlloc=zAlloc,nClusters=nClusters,hyperParams=hyperParams,
			wMat=wMat,nSubjects=nSubjects,
			nFixedEffects=nFixedEffects,
			nSweeps=nSweeps,nTableNames=nTableNames)
		pY <- -nSubjects * laplaceMode + 0.5*log(det(laplaceHessian)) -
			(nClusters+nFixedEffects)*0.5*log(nSubjects)+(nClusters+nFixedEffects)*0.5*log(2*pi)
		margModPost <- pY+pZpX  
		output<-list(margModPost=margModPost,laplacePar=laplaceOut$par[(nClusters+1):(nClusters+nFixedEffects)],pY=pY)
	} else {
		margModPost <- pZpX  
		output<-list(margModPost=margModPost)
	}
	return(output)
  
}

setHyperparams<-function(shapeAlpha=NULL,rateAlpha=NULL,aPhi=NULL,mu0=NULL,Tau0=NULL,R0=NULL,
	kappa0=NULL,nu0=NULL,muTheta=NULL,sigmaTheta=NULL,dofTheta=NULL,muBeta=NULL,sigmaBeta=NULL,dofBeta=NULL,
	shapeTauEpsilon=NULL,rateTauEpsilon=NULL,aRho=NULL,bRho=NULL,atomRho=NULL,shapeSigmaSqY=NULL,scaleSigmaSqY=NULL,
	rSlice=NULL,truncationEps=NULL,shapeTauCAR=NULL,rateTauCAR=NULL,shapeNu=NULL,scaleNu=NULL,initAlloc=NULL){
	out<-list()
	if (!is.null(shapeAlpha)){
		out$shapeAlpha<-shapeAlpha
	}
	if (!is.null(rateAlpha)){
		out$rateAlpha<-rateAlpha
	}
	if (!is.null(aPhi)){
		out$aPhi<-aPhi
	}
	if (!is.null(mu0)){
		out$mu0<-mu0
	}
	if (!is.null(Tau0)){
		out$Tau0<-Tau0
	}
	if (!is.null(R0)){
		out$R0<-R0
	}
	if (!is.null(kappa0)){
		out$kappa0<-kappa0
	}
	if (!is.null(nu0)){
		out$nu0<-nu0
	}
	if (!is.null(muTheta)){
		out$muTheta<-muTheta
	}
	if (!is.null(sigmaTheta)){
		out$sigmaTheta<-sigmaTheta
	}
	if (!is.null(dofTheta)){
		out$dofTheta<-dofTheta
	}
	if (!is.null(muBeta)){
		out$muBeta<-muBeta
	}
	if (!is.null(sigmaBeta)){
		out$sigmaBeta<-sigmaBeta
	}
	if (!is.null(dofBeta)){
		out$dofBeta<-dofBeta
	}
	if (!is.null(shapeTauEpsilon)){
		out$shapeTauEpsilon<-shapeTauEpsilon
	}
	if (!is.null(rateTauEpsilon)){
		out$rateTauEpsilon<-rateTauEpsilon
	}
	if (!is.null(aRho)){
		out$aRho<-aRho
	}
	if (!is.null(bRho)){
		out$bRho<-bRho
	}
	if (!is.null(atomRho)){
		out$atomRho<-atomRho
	}
	if (!is.null(shapeSigmaSqY)){
		out$shapeSigmaSqY<-shapeSigmaSqY
	}
	if (!is.null(scaleSigmaSqY)){
		out$scaleSigmaSqY<-scaleSigmaSqY
	}
	if (!is.null(rSlice)){
		out$rSlice<-rSlice
	}
	if (!is.null(truncationEps)){
		out$truncationEps<-truncationEps
	}
	if (!is.null(shapeTauCAR)){
	  out$shapeTauCAR<-shapeTauCAR
	}
	if (!is.null(rateTauCAR)){
	  out$rateTauCAR<-rateTauCAR
	}	
	if (!is.null(shapeNu)){
		out$shapeNu<-shapeNu
	}
	if (!is.null(scaleNu)){
		out$scaleNu<-scaleNu
	}
	if (!is.null(initAlloc)){
		out$initAlloc<-initAlloc
	}
	return(out)
}

	
# Compute Ratio of variances (for extra variation case)
computeRatioOfVariance<-function(runInfoObj){

	directoryPath=NULL
	extraYVar=NULL
	fileStem=NULL
	reportBurnIn=NULL
	nSweeps=NULL
	nFilter=NULL
	nSubjects=NULL
	nPredictSubjects=NULL
	nBurn=NULL
	
	for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])

	if (extraYVar==FALSE) stop("The ratio of variances can only be computed when extra variation in the response is included in the model.")

	# Construct the number of clusters file name
	nClustersFileName <- file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
	# Construct the allocation file name
	zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
	# Construct the allocation file name
	thetaFileName <- file.path(directoryPath,paste(fileStem,'_theta.txt',sep=''))
	# Construct the allocation file name
	epsilonFileName <- file.path(directoryPath,paste(fileStem,'_epsilon.txt',sep=''))
	
	# Restrict to sweeps after burn in
	firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
	lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter
	
	ratioOfVariance<-rep(0,length(lastLine-firstLine+1))
	for(sweep in firstLine:lastLine){
		currMaxNClusters<-scan(nClustersFileName,what=integer(),skip=sweep-1,n=1,quiet=T)
		zCurr<-1+scan(zFileName,what=integer(),skip=sweep-1,n=nSubjects+nPredictSubjects,quiet=T)
		zCurr<-zCurr[1:nSubjects]
		thetaCurr<-scan(thetaFileName,what=double(),skip=sweep-1,n=currMaxNClusters,quiet=T)
		thetaCurr<-thetaCurr[zCurr]
		vTheta<-var(thetaCurr)
		epsilonCurr<-scan(epsilonFileName,what=double(),skip=sweep-1,n=nSubjects,quiet=T)
		vEpsilon<-var(epsilonCurr)
		ratioOfVariance[sweep-firstLine+1]<-vTheta/(vTheta+vEpsilon)
		
	}
	return(ratioOfVariance)
	
}

# Function to convert a vector to a symmetric matrix (upper triangle)
vec2mat<-function (data = NA, nrow = 1) 
{
    nData <- length(data)

    nElem <- round(nrow * (nrow + 1)/2)

    result <- matrix(NA, nrow = nrow, ncol = nrow)

    result[lower.tri(result, diag = FALSE)] <- data
    result[upper.tri(result, diag = FALSE)] <- t(result)[upper.tri(result)]
    diag(result)<-0

    return(result)
}

# Function to plot a heatmap representing the dissimilarity matrix
# Some re-ordering of the observations is happening automatically
heatDissMat<-function(dissimObj,main=NULL,xlab=NULL,ylab=NULL)
{

	nSbj<-dissimObj$disSimRunInfoObj$nSubjects

	col.labels<-c("0","0.5","1")
	colours <- colorRampPalette(c("white","black"))(10)

	dissMat<-vec2mat(dissimObj$disSimMat,nrow=nSbj)
	heatmap(1-dissMat, keep.dendro=FALSE,symm=TRUE, Rowv=NA, labRow=FALSE, labCol=FALSE, margins=c(4.5,4.5), col= colours ,main = main, xlab=xlab, ylab=ylab)
	color.legend(0.95,0.7,1,1,legend = col.labels, colours, gradient="y",align="rb")

}

# function to plot the trace of the global parameters 
globalParsTrace<-function(runInfoObj, parameters = "nClusters",plotBurnIn=FALSE,whichBeta=1){

	directoryPath=NULL
	fileStem=NULL
	nSweeps=NULL
	nFilter=NULL
	nSubjects=NULL
	nBurn=NULL
	reportBurnIn=NULL
	alpha= NULL
	
	for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])

	if (parameters=="mpp") {
		parametersIn<-"margModPost"
	} else {
		parametersIn<-parameters
	}

	# Construct the parameter file name
	parFileName <- file.path(directoryPath,paste(fileStem,"_",parametersIn,".txt",sep=""))

	# read the data in
	parData<-read.table(parFileName)

	if(parameters== "nClusters") ylabPar<- "Number of clusters"
	if(parameters=="mpp") ylabPar<-"Log marginal model posterior"
	if(parameters=="beta") ylabPar<-"beta"
	if(parameters=="alpha") {
		if (alpha < -1) {
			ylabPar<-"alpha"
		} else {
			stop("ERROR: Parameter alpha is random only for the Dirichlet process prior with random alpha.") 
		}
	}

	xlabPars<-"Sweeps (after burn in)"


	if(plotBurnIn==TRUE){
		if (parameters=="mpp") stop("The mpp is only computed after the burn in. Set plotBurnIn=FALSE.")
		if (reportBurnIn==TRUE) {

			rangeSweeps<-1:((nBurn+nSweeps)/nFilter)
			rangeParData<-1:(length(rangeSweeps)/nFilter)
			xlabPars<-"Sweeps (including the burn in)"
		} else {
			stop("The function globalParsTrace cannot plot the burn in because the reportBurnIn option in profRegr was not set to TRUE and therefore the burn in has not been recorded.")
		}
	} else {
		if (parameters=="mpp"){
			firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
			lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter
			rangeSweeps<-firstLine:lastLine#(nBurn+1):(nBurn+nSweeps)
			rangeParData<-1:(nSweeps)	
		} else {
			if (reportBurnIn==TRUE) {
				firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
				lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter
				rangeSweeps<-firstLine:lastLine#(nBurn+1):(nBurn+nSweeps)
				rangeParData<-rangeSweeps
			} else {
				firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
				lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter
				rangeSweeps<-(nBurn/nFilter+1):((nBurn+nSweeps)/nFilter)
				rangeParData<-1:(nSweeps/nFilter)	

			}
		}
	}

	if (parameters=="alpha" || parameters=="nClusters"){
		plot(rangeSweeps,parData[rangeParData,1],type="l",lty=1,col="red",ylab=ylabPar,xlab=xlabPars,cex.axis=1.5,cex.lab=1.5)
	} else if (parameters=="mpp"){
		plot(1:nSweeps,parData[,1],type="l",lty=1,col="red",ylab=ylabPar,xlab=xlabPars,cex.axis=1.5,cex.lab=1.5)
	} else if (parameters=="beta"){
		plot(rangeSweeps,parData[rangeParData,whichBeta],type="l",lty=1,col="red",ylab=ylabPar,xlab=xlabPars,cex.axis=1.5,cex.lab=1.5)
	}

}

# plot posterior predictive densities
plotPredictions<-function(outfile,runInfoObj,predictions,logOR=FALSE){

	# ignores fixed effects
	
	nPredictedSubjects=NULL	
	directoryPath=NULL
	fileStem=NULL
	yModel=NULL
	xModel=NULL
	logOddsRatio=NULL
	weibullFixedShape=NULL

	for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])

	if (yModel!="Bernoulli"&&yModel!="Normal"&&yModel!="Survival") stop("This function has been developed for Bernoulli, Normal and Survival response only.")
	if (xModel=="Mixed") stop("This function has been developed for Discrete and Normal covariates only.")
	if (yModel=="Normal") logOR<-FALSE

	#if (runInfoObj$nFixedEffects>0) print("Note that fixed effects are not processed in this function.")
	
	predictResponseFileName = file.path(runInfoObj$directoryPath,paste(runInfoObj$fileStem,'_predict.txt',sep=''))
	relScenarios<-read.table(predictResponseFileName,header=FALSE,skip=1)

	# sistemare a seconda se logOR e' stato calcolato o no
	if (logOR==FALSE) {
		preds<-predictions$predictedYPerSweep[,,1]		
	} else {
		if (!is.null(predictions$logORPerSweep)){
			preds<-predictions$logORPerSweep
		} else {
			stop("Log OR (odds ratios) cannot be plotted because they have not been computed by calcPredictions. Re-run calcPredictions with option fullSweepLogOR=TRUE.")
		}
	}
	
	# output file
	pdf(outfile,onefile=TRUE)
   
	# Relevant scenarios
	nPredictSubjects<-runInfoObj$nPredictSubjects
	
	denObj<-vector(mode="list")
	for(i in 1:nPredictSubjects){
		denObj[[i]]<-density(na.omit(preds[,i]))#,bw=0.8) # removed over smoothing of predictions
	}
	
	for(k in 1:nPredictSubjects){
		plotDF<-data.frame('logOddsRatio'=denObj[[k]]$x,'density'=denObj[[k]]$y)
		plotObj<-ggplot(plotDF)
		plotObj<-plotObj+geom_line(aes(x=logOddsRatio,y=density),size=0.2)
		plotObj<-plotObj+theme(legend.position="none")
		plotObj<-plotObj+labs(x=ifelse(logOR==TRUE,"Log OR of response","Response"))+theme(axis.title.x=element_text(size=15))+labs(y="Density")+theme(axis.title.y=element_text(size=15,angle=90))
		plotObj<-plotObj+theme(axis.text.x=element_text(size=15))+theme(axis.text.y=element_text(size=15))
		print(plotObj)
	}
	output<-dev.off()
}

#.onAttach <- function(...) {
#	message("I need to collect case studies of the usage of the package PReMiuM. I would really appreciate if you could email me at liveranis@gmail.com if you use the package.")
# I need to collect case studies of the impact of my work (including developing and maintaining the package PReMiuM). I would really appreciate if you could email me at liveranis@gmail.com if you use the package.\n 
#	packageStartupMessage("Help support the development of PReMiuM: please share how PReMiuM is helping your work. Visit http://www.silvialiverani.com/support-premium/ or email liveranis@gmail.com.")
#\}

#.onLoad <- function(...) {
#	message("message from .onLoad via message")
#	packageStartupMessage("message from .onLoad via packageStartupMessage\n")
#}
