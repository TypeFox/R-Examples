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

# Generate simulated data for testing C++ PReMiuM
generateSampleDataFile<-function(clusterSummary){

	subjectsPerCluster<-clusterSummary$clusterSizes
	nSubjects<-sum(subjectsPerCluster)
	covariateType<-clusterSummary$covariateType
	nCovariates<-clusterSummary$nCovariates
	if (covariateType=="Mixed"){
		nDiscreteCovs<-clusterSummary$nDiscreteCovs
		nContinuousCovs<-clusterSummary$nContinuousCovs
	}
	missingDataProb<-clusterSummary$missingDataProb
	nFixedEffects<-clusterSummary$nFixedEffects
	nCategoriesY<-clusterSummary$nCategoriesY
	if (is.null(nCategoriesY)) nCategoriesY<-1
	includeCAR=clusterSummary$includeCAR
	nClusters=length(subjectsPerCluster)
   
	# Clustering covariates X
	X<-matrix(NA,nSubjects,nCovariates)

	k<-1

	#sample k for spatial CAR model
	if (includeCAR){
		VectoSample=c()
		for (k in 1:nClusters){
			VectoSample=c(VectoSample,rep(k,subjectsPerCluster[k]))
		}
		VectoSample=sample(VectoSample,nSubjects,replace=FALSE)
	}

	
	# Loop over subjects
	for(i in 1:nSubjects){
		if (includeCAR){
			clusterData<-clusterSummary$clusterData[[VectoSample[i]]]
		}else{
			if(i<=subjectsPerCluster[k]){
				clusterData<-clusterSummary$clusterData[[k]]
			}else{
				clusterData<-clusterSummary$clusterData[[k+1]]
				k<-k+1
				subjectsPerCluster[k]<-subjectsPerCluster[k]+subjectsPerCluster[k-1]
			}
		}
      
		# Loop over covariates to generate the X data
		if(covariateType=='Discrete'){
			for(j in 1:nCovariates){
				if(i>1&&runif(1)<missingDataProb){
					X[i,j]<-NA
				}else{
					u<-runif(1)
					nCategories<-length(clusterData$covariateProbs[[j]])
					for(kk in 1:nCategories){
						if(u<cumsum(clusterData$covariateProbs[[j]])[kk]){
							X[i,j]<-kk-1
							break
						}
					}
				}
			}
		}else if(covariateType=='Normal'){
			X[i,]<-clusterData$covariateMeans+t(chol(clusterData$covariateCovariance))%*%rnorm(nCovariates,0,1)
			for(j in 1:nCovariates){
				if(i>1&&runif(1)<missingDataProb){
					X[i,j]<-NA
				}
			}
		}else if(covariateType=='Mixed'){
			for(j in 1:nDiscreteCovs){
				if(i>1&&runif(1)<missingDataProb){
					X[i,j]<-NA
				}else{
					u<-runif(1)
					nCategories<-length(clusterData$covariateProbs[[j]])
					for(kk in 1:nCategories){
						if(u<cumsum(clusterData$covariateProbs[[j]])[kk]){
							X[i,j]<-kk-1
							break
						}
					}
				}
			}
			X[i,(nDiscreteCovs+1):nCovariates]<-clusterData$covariateMeans+t(chol(clusterData$covariateCovariance))%*%rnorm(nContinuousCovs,0,1)
			for(j in 1:nContinuousCovs){
				if(i>1&&runif(1)<missingDataProb){
					X[i,j]<-NA
				}
			}
		}
	}

	# Fixed Effects W
	if(nFixedEffects>0){
		W<-matrix(rnorm(nSubjects*nFixedEffects,0,1),nSubjects,nFixedEffects)
	}else{
		W<-NULL
	}

	# Spatial CAR term
	if (includeCAR) {
		tau=clusterSummary$TauCAR
		.write_neigh(nSubjects)
		M=.prec_Matrix()
		El=eigen(M)
		l=El$values
		e=El$vectors
		U=rnorm(nSubjects-1,0,1)
		U=U*(1/sqrt(l[1:(nSubjects-1)]))
		U=as.vector(e[,1:(nSubjects-1)]%*%matrix(U,ncol=1))
		U=U/sqrt(tau)    
	} else {
		U=rep(0,nSubjects)
	}

	# Response Vector Y
	Y<-rep(0,nSubjects)
	outcomeType<-clusterSummary$outcomeType
	if(nFixedEffects>0){
		if (outcomeType=='Categorical'){
			beta<-do.call(rbind,clusterSummary$fixedEffectsCoeffs)
		} else {
			beta<-as.matrix(clusterSummary$fixedEffectsCoeffs,nrow=1)
		}
	}
   	if (outcomeType=='Survival'){
		shape=clusterSummary$shape # shape parameter of the Weibull - constant across clusters
		censorT = clusterSummary$censorT # time to censor the outcomes
		event = clusterSummary$event # indicator variable for whether event has occured (1 = occured, 0 = censored)
	}

	if(outcomeType=='Poisson'){
		offset<-runif(nSubjects,clusterSummary$offsetLims[1],clusterSummary$offsetLims[2])
	}else{
		offset<-NULL
	}
   
	if(outcomeType=='Binomial'){	
		nTrials<-sample(clusterSummary$nTrialsLims[1]:clusterSummary$nTrials[2],nSubjects,replace=T)
	}else{
		nTrials<-NULL
	}

	if(outcomeType=='Normal'){
		sigmaSqY=clusterSummary$sigmaSqY
	}

	subjectsPerCluster<-clusterSummary$clusterSizes
	k<-1
	# Loop over subjects
	for(i in 1:nSubjects){
            
		if (includeCAR){
			theta<-clusterSummary$clusterData[[VectoSample[i]]]$theta
		}else{        
			if(i<=subjectsPerCluster[k]){
				theta<-clusterSummary$clusterData[[k]]$theta
				if (outcomeType=="Survival") shapeTmp<-clusterSummary$shape[k]
			}else{
				theta<-clusterSummary$clusterData[[k+1]]$theta
				if (outcomeType=="Survival") shapeTmp<-clusterSummary$shape[k+1]
				k<-k+1
				subjectsPerCluster[k]<-subjectsPerCluster[k]+subjectsPerCluster[k-1]
			}
		}
		mu<-theta
		if(nFixedEffects>0){
			if (outcomeType=='Categorical'){			
				for (kk in 2:nCategoriesY){
					mu[kk]<-mu[kk]+sum(beta[kk,]*W[i,])
				}
			} else {
				mu<-mu+sum(beta*W[i,])
			}
		} 
		if(outcomeType=='Poisson'){
			mu<-mu+U[i]
			mu<-mu+log(offset[i])
			Y[i]<-rpois(1,exp(mu))
		}else if(outcomeType=='Bernoulli'){
			p<-1/(1+exp(-mu))
			if(runif(1)<p){
				Y[i]<-1
			}else{
				Y[i]<-0
			}
		}else if(outcomeType=='Binomial'){
			p<-1/(1+exp(-mu))	
			Y[i]<-sum(runif(nTrials[i])<p)
		}else if(outcomeType=='Normal'){
			Y[i]<-rnorm(1,mu+U[i],sqrt(sigmaSqY))
		}else if (outcomeType=='Categorical'){
			p<-vector()
			sumMu<-sum(exp(mu))		
			p[1]<-1/sumMu
			for (kk in 2:nCategoriesY) p[kk]<-exp(mu[kk])/sumMu
			Y[i]<-which(rmultinom(1,1,p)==1)-1
		}else if (outcomeType == 'Survival'){
			Y[i] <-rWEI2(1, exp(mu), shapeTmp)#rlnorm(1,exp(mu),shapeTmp)#  #rweibull(1,exp(mu),shapeTmp)         #scale = exp(mu) 
			if (Y[i] >  censorT){  
				Y[i] <- censorT 
				event[i] <- 0
			} else {
				Y[i] <- Y[i]
				event[i] <- 1
			}
		}
	}

	# Write the output
	covNames<-paste('Variable',seq(1,nCovariates,1),sep="")
	if(nFixedEffects>0){
		fixEffNames<-paste('FixedEffects',seq(1,nFixedEffects,1),sep="")
	} else {
		fixEffNames=NULL
	}
	outData<-data.frame(cbind(matrix(Y),X))
	colnames(outData) <- c("outcome",covNames)
	out<-list(inputData=outData,covNames=covNames,xModel=covariateType,yModel=outcomeType)
	out$nCovariates <- nCovariates
	if (covariateType=="Mixed"){
		discreteCovs<-covNames[1:nDiscreteCovs]
		continuousCovs<-covNames[(nDiscreteCovs+1):nCovariates]
		out$nDiscreteCovs <- nDiscreteCovs
		out$nContinuousCovs <- nContinuousCovs
		out$discreteCovs <- discreteCovs
		out$continuousCovs <- continuousCovs
	}
	if(nFixedEffects>0){
		outData<-data.frame(cbind(outData,W))
		colnames(outData) <- c("outcome",covNames,fixEffNames)
		out$inputData <- outData
		out$fixedEffectNames <- fixEffNames
	}
	if(clusterSummary$outcomeType=="Poisson"){
		outData<-data.frame(cbind(outData,offset))
		out$inputData <- outData
		colnames(outData) <- c("outcome",covNames,fixEffNames,"outcomeT")
		out$inputData <- outData
		out$outcomeT <- "outcomeT"
	}
	if(clusterSummary$outcomeType=="Binomial"){
		outData<-data.frame(cbind(outData,nTrials))
		out$inputData <- outData	
		colnames(outData) <- c("outcome",covNames,fixEffNames,"outcomeT")
		out$inputData <- outData
		out$outcomeT <- "outcomeT"
	}
	if(clusterSummary$outcomeType=="Survival"){
		outData<-data.frame(cbind(outData,event))
		out$inputData <- outData	
		colnames(outData) <- c("outcome",covNames,fixEffNames, "event")
		out$inputData <- outData
		out$shape <- shape
		out$censorT <- censorT
	}
	if(clusterSummary$includeCAR){
		out$uCAR=U
		out$TauCAR=tau
		out$Permutation=VectoSample
	}
	return(out)
}


############################
# Sample datasets
clusSummaryCategoricalDiscrete<-function(){list(
	'outcomeType'='Categorical',
	'covariateType'='Discrete',
	'nCovariates'=5,
	'nCategories'=c(3,3,3,3,3),
	'nFixedEffects'=2,
	'fixedEffectsCoeffs'=list(c(0,0),c(-3,1),c(3,7)),
	'nCategoriesY'=3,
	'offsetLims'=c(0.9,1.1),
	'missingDataProb'=0.001,
	'nClusters'=3,
	'clusterSizes'=c(200,300,400),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=c(0,3,0.5),
		'covariateProbs'=list(c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1))),
		list('theta'=c(0,0.5,3),
		'covariateProbs'=list(c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.8,0.1,0.1))),
		list('theta'=c(0,-2,-2),
		'covariateProbs'=list(c(0.8,0.1,0.1),
			c(0.1,0.1,0.8),
			c(0.1,0.6,0.3),
			c(0.1,0.1,0.8),
			c(0.1,0.1,0.8)))))
}

clusSummaryBernoulliDiscrete<-function(){
	list(
	'outcomeType'='Bernoulli',
	'covariateType'='Discrete',
	'nCovariates'=5,
	'nCategories'=c(3,3,3,3,3),
	'nFixedEffects'=2,
	'fixedEffectsCoeffs'=c(0.1,-0.5),
	'missingDataProb'=0,
	'nClusters'=5,
	'clusterSizes'=c(200,200,200,200,200),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=log(9),
		'covariateProbs'=list(c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1))),
		list('theta'=log(2),
		'covariateProbs'=list(c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.1,0.8))),
		list('theta'=0,
		'covariateProbs'=list(c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1))),
		list('theta'=log(1/2),
		'covariateProbs'=list(c(0.1,0.1,0.8),
			c(0.1,0.8,0.1),
			c(0.8,0.1,0.1),
			c(0.1,0.1,0.8),
			c(0.8,0.1,0.1))),
		list('theta'=log(1/9),
		'covariateProbs'=list(c(0.1,0.1,0.8),
			c(0.1,0.1,0.8),
			c(0.1,0.1,0.8),
			c(0.1,0.1,0.8),
			c(0.1,0.1,0.8)))))
}

clusSummaryPoissonDiscrete<-function(){
	list(
	'outcomeType'='Poisson',
	'covariateType'='Discrete',
	'nCovariates'=5,
	'nCategories'=c(2,2,3,3,4),
	'nFixedEffects'=1,
	'fixedEffectsCoeffs'=c(0.01),
	'offsetLims'=c(0.9,1.1),
	'missingDataProb'=0.001,
	'nClusters'=5,
	'clusterSizes'=c(150,250,250,250,150),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=log(10),
		'covariateProbs'=list(c(0.8,0.2),
			c(0.2,0.8),
			c(0.6,0.2,0.2),
			c(0.25,0.5,0.25),
			rep(0.25,4))),
		list('theta'=log(5),
		'covariateProbs'=list(c(0.8,0.2),
			c(0.9,0.1),
			c(0.7,0.15,0.15),
			c(0.6,0.2,0.2),
			rep(0.25,4))),
		list('theta'=log(2.5),
		'covariateProbs'=list(c(0.2,0.8),
			c(0.9,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			rep(0.25,4))),
		list('theta'=log(1),
		'covariateProbs'=list(c(0.2,0.8),
			c(0.3,0.7),
			c(0.15,0.7,0.15),
			c(0.8,0.1,0.1),
			rep(0.25,4))),
		list('theta'=log(0.1),
		'covariateProbs'=list(c(0.2,0.8),
			c(0.5,0.5),
			c(0.1,0.1,0.8),
			c(0.8,0.1,0.1),
			c(0.1,0.1,0.1,0.7)))))
}

clusSummaryNormalDiscrete<-function(){
	list(
   'outcomeType'='Normal',
   'covariateType'='Discrete',
   'nCovariates'=5,
   'nCategories'=c(3,3,3,3,3),
   'nFixedEffects'=2,
   'fixedEffectsCoeffs'=c(0.1,-0.5),
   'sigmaSqY'=1,
   'missingDataProb'=0,
   'nClusters'=5,
   'clusterSizes'=c(500,700,300,600,900),
   'includeCAR'=FALSE,
   'TauCAR'=100,
   'clusterData'=list(list('theta'=10,
                           'covariateProbs'=list(c(0.8,0.1,0.1),
                                                 c(0.8,0.1,0.1),
                                                 c(0.8,0.1,0.1),
                                                 c(0.8,0.1,0.1),
                                                 c(0.8,0.1,0.1))),
                      list('theta'=5,
                           'covariateProbs'=list(c(0.8,0.1,0.1),
                                                 c(0.8,0.1,0.1),
                                                 c(0.1,0.8,0.1),
                                                 c(0.1,0.8,0.1),
                                                 c(0.1,0.1,0.8))),
                      list('theta'=0,
                           'covariateProbs'=list(c(0.1,0.8,0.1),
                                                 c(0.1,0.8,0.1),
                                                 c(0.1,0.8,0.1),
                                                 c(0.1,0.8,0.1),
                                                 c(0.1,0.8,0.1))),
                      list('theta'=-5,
                           'covariateProbs'=list(c(0.1,0.1,0.8),
                                                 c(0.1,0.8,0.1),
                                                 c(0.8,0.1,0.1),
                                                 c(0.1,0.1,0.8),
                                                 c(0.8,0.1,0.1))),
                      list('theta'=-10,
                           'covariateProbs'=list(c(0.1,0.1,0.8),
                                                 c(0.1,0.1,0.8),
                                                 c(0.1,0.1,0.8),
                                                 c(0.1,0.1,0.8),
                                                 c(0.1,0.1,0.8)))))
}

clusSummaryPoissonNormal<-function(){
	list(
	'outcomeType'='Poisson',
	'covariateType'='Normal',
	'nCovariates'=2,
	'nFixedEffects'=2,
	'fixedEffectsCoeffs'=c(-0.05,0.1),
	'offsetLims'=c(0.9,1.1),
	'missingDataProb'=0.001,
	'nClusters'=3,
	'clusterSizes'=c(500,600,400),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=log(10),
			'covariateMeans'=c(0,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
		list('theta'=log(3),
			'covariateMeans'=c(3,2),
			'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
		list('theta'=log(0.1),
			'covariateMeans'=c(10,-5),
			'covariateCovariance'=matrix(c(2,0.7,0.7,1),nrow=2))))
}

clusSummaryPoissonNormalSpatial<-function(){
  list(
    'outcomeType'='Poisson',
    'covariateType'='Normal',
    'nCovariates'=2,
    'nFixedEffects'=2,
    'fixedEffectsCoeffs'=c(-0.05,0.1),
    'offsetLims'=c(0.9,1.1),
    'missingDataProb'=0.001,
    'nClusters'=3,
    'clusterSizes'=c(50,60,40),
    'includeCAR'=TRUE,
    'TauCAR'=100,
    'clusterData'=list(list('theta'=log(10),
                            'covariateMeans'=c(0,2),
                            'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
                       list('theta'=log(3),
                            'covariateMeans'=c(3,2),
                            'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
                       list('theta'=log(0.1),
                            'covariateMeans'=c(10,-5),
                            'covariateCovariance'=matrix(c(2,0.7,0.7,1),nrow=2))))
}


clusSummaryNormalNormalSpatial<-function(){
  list(
    'outcomeType'='Normal',
    'covariateType'='Normal',
    'nCovariates'=2,
	'sigmaSqY'=1,
    'nFixedEffects'=2,
    'fixedEffectsCoeffs'=c(-0.05,0.1),
    'missingDataProb'=0.001,
    'nClusters'=3,
    'clusterSizes'=c(50,60,40),
    'includeCAR'=TRUE,
    'TauCAR'=100,
    'clusterData'=list(list('theta'=log(10),
                            'covariateMeans'=c(0,2),
                            'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
                       list('theta'=log(3),
                            'covariateMeans'=c(3,2),
                            'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
                       list('theta'=log(0.1),
                            'covariateMeans'=c(10,-5),
                            'covariateCovariance'=matrix(c(2,0.7,0.7,1),nrow=2))))
}


clusSummaryBinomialNormal<-function(){
	list(
	'outcomeType'='Binomial',
	'covariateType'='Normal',
	'nCovariates'=2,
	'nFixedEffects'=1,
	'fixedEffectsCoeffs'=c(0.1),
	'nTrialsLims'=c(5,15),
	'missingDataProb'=0.001,
	'nClusters'=3,
	'clusterSizes'=c(200,700,100),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=log(10),
			'covariateMeans'=c(0,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
		list('theta'=log(3),
			'covariateMeans'=c(3,2),
			'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
		list('theta'=log(0.1),
			'covariateMeans'=c(10,-5),
			'covariateCovariance'=matrix(c(2,0.7,0.7,1),nrow=2))))
}
	
clusSummaryNormalNormal<-function(){
	list(
	'outcomeType'='Normal',
	'covariateType'='Normal',
	'nCovariates'=2,
	'nFixedEffects'=0,
	'sigmaSqY'=1,
	'missingDataProb'=0,
	'nClusters'=3,
	'clusterSizes'=c(300,500,400),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=-5,
			'covariateMeans'=c(0,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
		list('theta'=0,
			'covariateMeans'=c(3,2),
			'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
		list('theta'=5,
			'covariateMeans'=c(10,-5),
			'covariateCovariance'=matrix(c(2,0.9,0.9,1),nrow=2))))
}

clusSummaryVarSelectBernoulliDiscrete<-function(){
	list(
	'outcomeType'='Bernoulli',
	'covariateType'='Discrete',
	'nCovariates'=10,
	'nCategories'=rep(2,10),
	'nFixedEffects'=0,
	'missingDataProb'=0,
	'nClusters'=5,
	'clusterSizes'=c(200,200,200,200,200),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=log(1.0/9.0),
	'covariateProbs'=list(c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.9,0.1),
			c(0.1,0.9))),
		list('theta'=log(3.0/7.0),
		'covariateProbs'=list(c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.1,0.9))),
		list('theta'=0,
		'covariateProbs'=list(c(0.9,0.1),
			c(0.1,0.9),
			c(0.9,0.1),
			c(0.1,0.9),
			c(0.9,0.1),
			c(0.1,0.9),
			c(0.9,0.1),
			c(0.1,0.9),
			c(0.9,0.1),
			c(0.1,0.9))),
		list('theta'=log(7.0/3.0),
		'covariateProbs'=list(c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.9,0.1),
			c(0.1,0.9))),
		list('theta'=log(9.0),
		'covariateProbs'=list(c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.1,0.9),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.9,0.1),
			c(0.1,0.9)))))
}

clusSummaryBernoulliMixed<-function(){
	list(
	'outcomeType'='Bernoulli',
	'covariateType'='Mixed',
	'nCovariates'=5,
	'nDiscreteCovs'=3,
	'nContinuousCovs'=2,
	'nCategories'=c(3,3,3),
	'nFixedEffects'=0,
	'missingDataProb'=0,
	'nClusters'=3,
	'clusterSizes'=c(300,300,300),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=log(10),
		'covariateProbs'=list(c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1)),
			'covariateMeans'=c(-10,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
		list('theta'=log(3),
		'covariateProbs'=list(c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1)),
			'covariateMeans'=c(3,20),
			'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
		list('theta'=log(0.1),
		'covariateProbs'=list(c(0.1,0.1,0.8),
			c(0.1,0.1,0.8),
			c(0.1,0.1,0.8)),
			'covariateMeans'=c(10,-5),
			'covariateCovariance'=matrix(c(2,0.9,0.9,1),nrow=2))))
}

			
clusSummaryBernoulliDiscreteSmall<-function(){
	list(
	'outcomeType'='Bernoulli',
	'covariateType'='Discrete',
	'nCovariates'=5,
	'nCategories'=c(3,3,3,3,3),
	'nFixedEffects'=2,
	'fixedEffectsCoeffs'=c(0.1,-0.5),
	'missingDataProb'=0,
	'nClusters'=5,
	'clusterSizes'=c(100,100,100,100,100),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=log(9),
		'covariateProbs'=list(c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1))),
		list('theta'=log(2),
		'covariateProbs'=list(c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.1,0.8))),
		list('theta'=0,
		'covariateProbs'=list(c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1),
			c(0.1,0.8,0.1))),
		list('theta'=log(1/2),
		'covariateProbs'=list(c(0.1,0.1,0.8),
			c(0.1,0.8,0.1),
			c(0.8,0.1,0.1),
			c(0.1,0.1,0.8),
			c(0.8,0.1,0.1))),
		list('theta'=log(1/9),
		'covariateProbs'=list(c(0.1,0.1,0.8),
			c(0.1,0.1,0.8),
			c(0.1,0.1,0.8),
			c(0.1,0.1,0.8),
			c(0.1,0.1,0.8)))))
}

			
clusSummaryBernoulliNormal<-function(){
	list(
	'outcomeType'='Bernoulli',
	'covariateType'='Normal',
	'nCovariates'=2,
	'nFixedEffects'=0,
	'missingDataProb'=0,
	'nClusters'=5,
	'clusterSizes'=c(100,100,100,100,100),
	'includeCAR'=FALSE,
	'TauCAR'=100,
	'clusterData'=list(list('theta'=log(9),
			'covariateMeans'=c(0,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
		list('theta'=log(2),
			'covariateMeans'=c(0,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
		list('theta'=0,
			'covariateMeans'=c(0,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
		list('theta'=log(1/2),
			'covariateMeans'=c(0,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),
		list('theta'=log(1/9),
			'covariateMeans'=c(0,2),
			'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2))))

}


clusSummaryWeibullDiscrete<-function(){                 
	list(
	'outcomeType'='Survival',
	'covariateType'='Discrete',
	'nCovariates'=5,
	'nCategories'=c(2,2,3,3,4),
	'nFixedEffects'=1,                                      
	'fixedEffectsCoeffs'=c(0),
	'shape'=c(2,2.4,3),
	'censorT'=50,                                                                  
	'missingDataProb'=0,
	'nClusters'=3,                                      
	'clusterSizes'=c(250,250,250),
	'includeCAR'=FALSE,
	'clusterData'=list(list('theta'=log(0.01),                   
		'covariateProbs'=list(c(0.8,0.2),
			c(0.8,0.2),
			c(0.8,0.1,0.1),
			c(0.8,0.1,0.1),
			rep(0.25,4))),
		list('theta'=log(0.001),            
		'covariateProbs'=list(c(0.8,0.2),
			c(0.2,0.8),
			c(0.1,0.1,0.8),
			c(0.1,0.8,0.1),
			rep(0.25,4))),
		list('theta'= log(0.005),  
		'covariateProbs'=list(c(0.2,0.8),
			c(0.2,8),
			c(0.1,0.1,0.8),
			c(0.1,0.8,0.1),
			c(0.1,0.1,0.1,0.7)))))

}








############# Additional functions for spatial term


#Function for mapping the spatial effect when data are generated
mapforGeneratedData=function(u, del=NULL, palette='RGB', main=''){
	#Fill the matrix with u
	if (floor(sqrt(length(u)))==sqrt(length(u))){n_u=length(u)
	} else {n_u=(floor(sqrt(length(u)))+1)^2}
	u_compl=rep(NA, n_u)
	u_compl[1:length(u)]=u
	n=sqrt(n_u)
	mat=matrix(NA,nrow=n,ncol=n)
	for (i in 1:n){mat[i,]=u[((i-1)*n)+1:n]}
	lagS=seq(0:n)
	lagT=seq(0:n)

	#Set plot parameters
	if (palette=='BW'){paletteBleue<-colorRampPalette(c("Gray 92","Gray 16"))
	}else{paletteBleue<-colorRampPalette(c("blue","cyan","yellow","orange","red"))}
	if(is.null(del)){
		deltot=seq(range(u,na.rm=TRUE)[1],range(u,na.rm=TRUE)[2],length.out=100)
		deltot[1]<-deltot[1]-(deltot[2]-deltot[1])/100
		deltot[length(deltot)]<-deltot[length(deltot)]+(deltot[length(deltot)]-deltot[length(deltot)-1])/100  
		labs=c(seq(1,length(deltot),20),length(deltot))
	}else{
		labs=c(1)
		deltot=c()
		for(i in 1:(length(del)-1)){
			deltot=c(deltot,seq(del[i],del[i+1],length.out=20))
			if(deltot[length(deltot)]==del[i+1]){deltot=deltot[-length(deltot)]}
			labs=c(labs,length(deltot)+1)
		}
		deltot=c(deltot,del[length(del)])
	}
	pardef=par(no.readonly=TRUE)
	par(mar=c(4,4,4,8),tck=0.02,mgp=c(3,0.2,0),las=1)
       
	#Begin plot      
	plot.new()
	plot.window(xlim=range(lagS),ylim=range(lagT),xlab="oui", ylab="non")
	for ( s in 1:length(mat[1,])){
		for ( t in 1:length(mat[,1])){
			rect(xleft=lagS[s],xright=lagS[s+1],ybottom=lagT[t],ytop=lagT[t+1],col=paletteBleue(length(deltot-1))[findInterval(mat[t,s],deltot,all.inside=TRUE)],border=NA)
    		}
	}
                    
	axis(1, at=lagS[seq(1,length(lagS),2)],labels=as.character(lagS[seq(1,length(lagS),2)]),tick=TRUE, pos=c(lagS[1],lagT[1]), padj=0, cex.axis=1)
	axis(2, at=lagT[seq(1,length(lagT),2)],labels=as.character(lagT[seq(1,length(lagT),2)]),tick=TRUE,pos=c(lagS[1],lagT[1]), cex.axis=1)
	axis(3, at=lagS[seq(1,length(lagS),2)],labels=FALSE,tick=TRUE, pos=c(lagT[length(lagT)],lagS[1]))
	axis(4,at=lagT[seq(1,length(lagT),2)],labels=FALSE,tick=TRUE, pos=c(lagS[length(lagS)],lagT[length(lagT)]))
	rect(xleft=range(lagS)[1],xright=range(lagS)[2],ybottom=range(lagT)[1],ytop=range(lagT)[2],border='black',lwd=1)        
	mtext(main, side=3, line=0.5, cex=2)
  
	colors=.color.scale(paletteBleue(length(deltot)-1),deltot,name="", unit="",labels=labs)
	.cs.draw(colors,border=NA,horiz=F, cex=1,digits=2, side=1, length=0.8, offset=0, pos=1, width = 0.06)
	par(pardef)
}


########################################################################################################
#function for producing text file with lattice neighbourhood structure that may be read 

.write_neigh=function(nSubjects, file=NULL){
  if (is.null(file)) {file='Neighbours.txt'}
  if (floor(sqrt(nSubjects))==sqrt(nSubjects)){ n=sqrt(nSubjects)
  }else{ n=floor(sqrt(nSubjects))+1 }
  out=c(nSubjects)
  for ( i in 1:nSubjects){
    if (i <=n){
      if (i==1) out=c(out, paste(i,2,2,1+n))
      else if (i==n) out=c(out, paste(i,2,i-1,2*n))
      else out=c(out, paste(i, 3, i-1,i+1,i+n))
    } else if (i <=(nSubjects-n)){
      if (floor((i-1)/n)==(i-1)/n) out=c(out, paste(i,3,i-n,i+1,i+n))
      else if (floor(i/n)==i/n) out =c(out, paste(i,3,i-n,i-1,i+n))
      else out =c(out, paste (i,4,i-n,i-1,i+1,i+n))
    } else {
      if (floor((i-1)/n)==(i-1)/n) out=c(out, paste(i,2,i-n,i+1))
      else if (floor(i/n)==i/n) out =c(out, paste(i,2,i-n,i-1))
      else if (i==nSubjects) {
        if (floor((i-1)/n)==(i-1)/n) out=c(out,paste(i,1,i-n))
        else out=c(out, paste(i,2,i-n,i-1))
      } else {out =c(out, paste(i,3,i-n,i-1,i+1))}      
    }
  }
  writeLines(out,file)
}

.prec_Matrix=function(file=NULL){
  if (is.null(file)) file='Neighbours.txt'
  con=file(file)
  open(con)
  nSubjects=as.numeric(readLines(con,n=1,warn=FALSE))
  MAT=matrix(0,ncol=nSubjects,nrow=nSubjects)
  for (i in 1:nSubjects){
    tmp=as.numeric(strsplit(readLines(con,n=1,warn=FALSE),split=' ')[[1]])
    MAT[i,i]=tmp[2]
    MAT[i,tmp[2+1:tmp[2]]]=-1
  }
  close(con)
  return(MAT)
}

.cs.draw <- function (color.scale, name = NULL, unit = NULL, length = 0.8,
    width = 0.03, horiz = T, pos = 1.09, side = if (pos > 0.5) -1 else 1,
    cex = 1, offset = 0, border = NULL, lty = NULL, lwd = par("lwd"),
    xpd = T, digits = 4, roundfunc = zapsmall)
{
    nc <- length(color.scale$cols)
    if (is.null(name))
        name <- color.scale$name
    if (is.null(unit))
        unit <- color.scale$unit
    if (horiz) {
        xc <- (0.5 - offset) * par("usr")[1] + (0.5 + offset) *
            par("usr")[2]
        xd <- length * (par("usr")[2] - par("usr")[1])
        x1 <- xc - xd/2
        x2 <- xc + xd/2
        x <- seq(x1, x2, , nc + 1)
        ya <- par("usr")[4] - par("usr")[3]
        yd <- width * (par("usr")[4] - par("usr")[3])
        y1 <- par("usr")[3] + pos * ya
        y2 <- y1 + side * yd
        for (i in 1:nc) {
            rect(x[i], y1, x[i + 1], y2, col = color.scale$cols[i],
                border = border, lty = lty, lwd = lwd, xpd = xpd)
            if (i %in% color.scale$labels)
                text(x[i], y2, adj = c(0.5, -0.75 * side + 0.5),
                  roundfunc(color.scale$breaks[i], digits = digits),
                  xpd = xpd, cex = cex)
        }
        text(x[nc + 1], y2, adj = c(0.5, -0.75 * side + 0.5),
            roundfunc(color.scale$breaks[nc + 1], digits = digits),
            xpd = xpd, cex = cex)
        text(xc, y1, adj = c(0.5, 0.75 * side + 0.5), name, xpd = xpd,
            cex = cex)
        text(x[nc + 1], (y1 + y2)/2, adj = c(-0.05, 0.5), unit,
            xpd = xpd, cex = cex)
    }
    else {
        yc <- (0.5 - offset) * par("usr")[3] + (0.5 + offset) *
            par("usr")[4]
        yd <- length * (par("usr")[4] - par("usr")[3])
        y1 <- yc - yd/2
        y2 <- yc + yd/2
        y <- seq(y1, y2, , nc + 1)
        xa <- par("usr")[2] - par("usr")[1]
        xd <- width * (par("usr")[2] - par("usr")[1])
        x1 <- par("usr")[1] + pos * xa
        x2 <- x1 + side * xd
        for (i in 1:nc) {
            rect(x1, y[i], x2, y[i + 1], col = color.scale$cols[i],
                border = border, lty = lty, lwd = lwd, xpd = xpd)
            if (i %in% color.scale$labels)
                text(x2, y[i], adj = c(-0.75 * side + 0.5, 0.5),
                  roundfunc(color.scale$breaks[i], digits = digits),
                  xpd = xpd, cex = cex)
        }
        text(x2, y[nc + 1], adj = c(-0.75 * side + 0.5, 0.5),
            roundfunc(color.scale$breaks[nc + 1], digits = digits),
            xpd = xpd, cex = cex)
        text(x1, yc, adj = c(0.5, -0.75 * side + 0.5), name,
            xpd = xpd, srt = 90, cex = cex)
        text(x1, y[1], adj = c(-0.75 * side + 0.5, 1.5), unit,
            xpd = xpd, cex = cex)
    }
} 

#Function to create color.scale
.color.scale=function(vectofcols, vectofbreaks, name="name", unit="", labels){
  res=list()
  res$cols=vectofcols
  res$breaks=vectofbreaks
  res$name=name
  res$unit=unit
  if (is.null(labels)){res$labels=seq(1,length(vectofbreaks),1)
  }else{res$labels=labels}
  return(res)
}


