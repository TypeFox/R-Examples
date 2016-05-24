

.womActiveSaveTab<-function(status=TRUE)
{
	theWidget("savePageBox")$setSensitive(status)	
}


.save_clicked_change_working_directory <- function(button)
{
	fileChoose(action="setwd", type="selectdir", text="Select a directory...")
	workingDir = getwd()
	theWidget("saveWorkingDirectoryEntry")$setText(workingDir)
}


.save_clicked_save_button <- function(button)
{
	analysisName = .womFileSafeString(theWidget("saveAnalysisNameEntry")$getText())	
	.womSaveAnalysis(analysisName)
}

.womSaveAnalysis <- function(analysisName)
{
	freezeGUI()
	on.exit(thawGUI())
	filename = paste(analysisName,".Rdata",sep="")
	loaded = try(save(file=filename,envir=wommbatAnalysis), silent=T)
	if (inherits(loaded, "try-error")){
		gWidgets::gmessage(paste("The R data file could not be saved. Is the location read-only?.",
						sep=""), title="File error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)
	}

	if(theWidget("saveTextFilesCheckbox")$getActive())
	{
		.womExportAnalysis(analysisName)
		.womSetStatusBarText("Saved R file and exported text files.")
	}else
	{
		if(length(wommbatAnalysis$Models))
		{
			.womSetStatusBarText("Saved R file.")
		}
	}
	
}


.womExportAnalysis <- function(analysisName="analysis")
{
	directoryName = paste(analysisName,"_textfiles",sep="")
	if(file.exists(directoryName)) unlink(directoryName)
	success = dir.create(directoryName)
	if (!success){
			gWidgets::gmessage(paste("The analysis text file directory could not be created. Is the location read-only?.",
						sep=""), title="File error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(0)
		}
	.womCreateSummaryFile(analysisName,directoryName)
	.womCreateModelFiles(analysisName,directoryName)

	write.csv(file=file.path(directoryName,paste(analysisName,"_data.csv",sep="")),wommbatAnalysis$data)
}

.womCreateModelFiles <- function(analysisName,directoryName)
{
		models = names(wommbatAnalysis$Models)
		for(modelName in models)
		{
			subDirectoryName = .womFileSafeString(modelName)
			subPath = file.path(directoryName,subDirectoryName)
			success = dir.create(subPath)
			if (!success){
				gWidgets::gmessage(paste("One of the analysis text file subdirectories could not be saved. Is the location read-only?.",
						sep=""), title="File error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
				return(0)
			}
			.womCreateModelSummaryFile(modelName,subPath)
			write.csv(file=file.path(subPath,paste(.womFileSafeString(modelName),"_ContinuousData.csv",sep="")),wommbatAnalysis$Models[[modelName]]$model$newDat2Cont)
			write.csv(file=file.path(subPath,paste(.womFileSafeString(modelName),"_CategoricalData.csv",sep="")),wommbatAnalysis$Models[[modelName]]$model$newDat2Cat)
			write.csv(file=file.path(subPath,paste(.womFileSafeString(modelName),"_pointEstimates.csv",sep="")),wommbatAnalysis$Models[[modelName]]$results$pointEst)
			.womCreateModelChainFiles(modelName,subPath)
		}
}

.womCreateModelChainFiles <- function(modelName,subPath)
{
	effectChainsFilename = paste(.womFileSafeString(modelName),"_effectChains.csv",sep="")
	covChainsFilename = paste(.womFileSafeString(modelName),"_covChains.csv",sep="")
	corChainsFilename = paste(.womFileSafeString(modelName),"_corChains.csv",sep="")
	meanChainsFilename = paste(.womFileSafeString(modelName),"_meanChains.csv",sep="")
	myModel = wommbatAnalysis$Models[[modelName]]
	
	if(is.null(myModel$results)) return(0)
	
	# write effect chains
	write.csv(file=file.path(subPath,effectChainsFilename),myModel$results$effectChains)
	
	# write covariance/correlation chains
	if(myModel$model$covNgroups>0)
	{
		covChains = corChains = NULL
		for(i in 1:myModel$model$covNgroups)
		{
			covChainsTemp = myModel$results$covChains[[i]]
			covChainsLabs = outer(dimnames(covChainsTemp)[[1]],dimnames(covChainsTemp)[[2]],paste,sep=" and ")
			covKeep =  which(lower.tri(covChainsTemp[,,1],diag=TRUE))
			corChainsTemp = myModel$results$covChains[[i]]
			corKeep =  which(lower.tri(corChainsTemp[,,1],diag=FALSE))
			dim(covChainsTemp) = c(myModel$model$covSizes[i]^2,myModel$settings$effectiveIters)
			dim(corChainsTemp) = dim(covChainsTemp)
			rownames(corChainsTemp) = covChainsLabs
			rownames(covChainsTemp) = covChainsLabs
			covChains = cbind(covChains,t(covChainsTemp)[,covKeep])
			corChains = cbind(covChains,t(corChainsTemp)[,corKeep])
		}
		write.csv(file=file.path(subPath,covChainsFilename),covChains)
		write.csv(file=file.path(subPath,corChainsFilename),corChains)
	}
	# write mean chains?
	
}

.womCreateModelSummaryFile<-function(modelName,subPath)
{
	filename = paste(.womFileSafeString(modelName),".txt",sep="")
	myModel = wommbatAnalysis$Models[[modelName]]
	
	if(wommbatAnalysis$Ktype==0){
		capMod = "Cowan"
	}else if(wommbatAnalysis$Ktype==1){
		capMod = "Pashler"
	}else if(wommbatAnalysis$Ktype==2){
		capMod = "Morey (Experimental)"
	}
	
	useCovMod = ifelse(myModel$model$covNgroups>0,"Yes","No")
	
	columns = paste(paste(wommbatAnalysis$SelCols[,1],wommbatAnalysis$SelCols[,2],sep=": "),collapse="\n")
	
	onK = paste("kappa    = ",paste(c("muK",as.character(myModel$model$SelEffs[as.logical(myModel$model$SelEffs[,5]),4])),collapse=" + "),collapse="")
	onZ = paste("logit(Z) = ",paste(c("muZ",as.character(myModel$model$SelEffs[as.logical(myModel$model$SelEffs[,6]),4])),collapse=" + "),collapse="")
	onG = paste("logit(G) = ",paste(c("muG",as.character(myModel$model$SelEffs[as.logical(myModel$model$SelEffs[,7]),4])),collapse=" + "),collapse="")

	if(myModel$model$covNgroups>0)
	{
		covs = vector(length=myModel$model$covNgroups)

		for(i in 1:myModel$model$covNgroups)
		{
			covs[i] = paste(dimnames(myModel$results$covChains[[i]])[[1]],sep=", ")
		}
		covs=paste("  Covariance Matrices:\n",paste(covs,collapse="\n"))
		wishartdf = paste("Wishart df     : ",myModel$priors$invWishartScalar,"\n",sep="")

	}else{
		covs=""
		wishartdf = ""
	}
	
	if(!is.null(myModel$results)){
		DICtext = paste(round(myModel$results$DIC,1)," (Effective parameters: ",round(myModel$results$pD,1),")",sep="")
		if(myModel$settings$useMH)
		{
			MCMCtext = paste(
				"\nMCMC Settings\n-------------------\n",
				"Type           : Random Walk Metropolis-Hastings\n",
				"Iterations     : ",myModel$settings$MCMCIters,"\n",
				"Burnin         : ",myModel$settings$burninIters,"\n",
				"Thin           : ",myModel$settings$MHThin,"\n",
				"Effective Iters: ",myModel$settings$effectiveIters,"\n",
				"Scale          : ",myModel$settings$MHscale,"\n\n",
			sep="")
		}else
		{
			MCMCtext = paste(
				"\nMCMC Settings\n-------------------\n",
				"Type           : Hybrid\n",
				"Iterations     : ",myModel$settings$MCMCIters,"\n",
				"Burnin         : ",myModel$settings$burninIters,"\n",
				"Epsilon        : (",myModel$settings$epsLow,", ",myModel$settings$epsUpp,")\n",
				"Leapfrog steps : ",myModel$settings$leapfrog,"\n\n",
			sep="")
		}
	}else{
		MCMCtext = "Model has not been analyzed yet.\n\n"
		DICtext = "NA\n"
	}
	  ## Write info
	outputInfo = paste(
		"WoMMBAT Analysis\n",
		"--------------------\n",
		"Analysis name: ",modelName,"\n",
		"Analysis time: ",as.character(myModel$results$timeAnalyzed),"\n",
		"Model Type   : ",capMod,"\n",
		"Covariances? : ",useCovMod,"\n",
		"\nPrior Specification\n-------------------\n",
		"Inverse Gamma  : (a=",myModel$priors$IGa0,", b=",myModel$priors$IGb0,")\n",
		"muK            : Normal(mu=",myModel$priors$muKmean,", sigma=",myModel$priors$muKSD,")\n",
		"muZ            : Normal(mu=",myModel$priors$muZmean,", sigma=",myModel$priors$muZSD,")\n",
		"muG            : Normal(mu=",myModel$priors$muGmean,", sigma=",myModel$priors$muGSD,")\n",
		wishartdf, MCMCtext,

		"\nModel\n-------------------\n",
		columns,"\n\n",
		onK,"\n",onZ,"\n",onG,"\n",
		covs,"\n\n",
		"DIC            : ",DICtext,"\n",
    sep="")


	touch = try(cat(outputInfo,file=file.path(subPath,filename),append=FALSE), silent=T)
	if (inherits(touch, "try-error")){
		gWidgets::gmessage(paste("The model summary file for '", modelName,"' could not be saved. Is the location read-only?.",
					sep=""), title="File error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)
	}
}

.womCreateSummaryFile <- function(analysisName,directoryName)
{
		header = paste("Model name",
					   "Has results?",
					   "Time analyzed",
					   "DIC",
					   "pD",
					   "Effective MCMC iterations",
					"\n",sep="\t")
		filename = file.path(directoryName,paste(analysisName,".txt",sep=""))
		touch = try(cat(header,file=filename,append=FALSE), silent=T)
		if (inherits(touch, "try-error")){
			gWidgets::gmessage(paste("The analysis summary file could not be saved. Is the location read-only?.",
						sep=""), title="File error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(0)
		}
				
		models = names(wommbatAnalysis$Models)
		for(modelName in models)
		{
			myModel = wommbatAnalysis$Models[[modelName]]
			cat(modelName,
				!is.null(myModel$results),
				as.character(myModel$results$timeAnalyzed),
				myModel$results$DIC,
				myModel$results$pD,
				myModel$settings$effectiveIters,
				"\n",file=filename,sep="\t",append=TRUE)
		}
}


