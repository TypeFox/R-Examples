##########################################################################
#gsAnalysis
#
# method to define a "gsAnalysis" 
##########################################################################
gsAnalysis <- function(
		name,
		gls = NULL,
		glsParameterNames = NULL,
		transformation = NULL,
		transformationParameterNames = NULL,
		gss = NULL,
		gssParameterNames = NULL,
		globalStat = NULL,
		globalStatParameterNames = NULL,
		significance = NULL,
		significanceParameterNames = NULL,
		testAlternative = c("greater", "less")){

	analysis <- list(
		name							= name,
		gls								= gls,
		glsParameterNames 				= glsParameterNames,
		transformation 					= transformation,
		transformationParameterNames 	= transformationParameterNames,
		gss 							= gss,
		gssParameterNames				= gssParameterNames,
		globalStat 						= globalStat,
		globalStatParameterNames		= globalStatParameterNames,
		significance 					= significance,
		significanceParameterNames		= significanceParameterNames,
		testAlternative					= testAlternative)

	class(analysis) <- "gsAnalysis"
	return(analysis)
}

##########################################################################
#print.gsAnalysis
##########################################################################
print.gsAnalysis <- function(x, ...){
	cat(x$name, "\n")
	if(is.null(x$globalStat)){
		cat("Gene level statistic:    ", x$gls,
			"\n                         (parameters: ",
		  	paste(x$glsParameterNames,collapse=","),")\n", sep="")
		cat("Transformation:          ", x$transformation,
			"\n                         (parameters: ",
			paste(x$transformationParameterNames,collapse=","),")\n", sep="")
		cat("Gene set statistic:      ", x$gss,
			"\n                         (parameters: ",
			paste(x$gssParameterNames,collapse=","),")\n", sep="")
		cat("Significance assessment: ", x$significance,
			"\n                         (parameters: ",
			paste(x$significanceParameterNames,collapse=","),")\n", sep="")
		cat("Test alternative:        ", x$testAlternative, "\n", sep ="")
	}else{
		cat("Global analysis:         ", x$globalStat,
			"\n                         (parameters: ",
			paste(x$globalStatParameterNames,collapse=","),")\n", sep="")
	}
}

##########################################################################
# Gene Set Enrichment Analysis: Subramanian et al. 2005
##########################################################################
analysis.gsea <- function(){
	return(gsAnalysis(
		name 							= "gsea",
		gls 							= "gls.cor",
		glsParameterNames 				= c("labs","method"),
		transformation 					= "transformation.abs",
		transformationParameterNames 	= NULL,
		gss 							= "gss.enrichmentScore",
		gssParameterNames 				= c("p"),
		globalStat 						= NULL,
		globalStatParameterNames 		= NULL,
		significance 					= "significance.permutation",
		significanceParameterNames 		= c("numSamples","labs"),
		testAlternative					= "greater"))
}

##########################################################################
# Fisher exact test for 2 sets of genes (coreset and geneset)
##########################################################################
analysis.customOverrepresentation <- function(){
	return(gsAnalysis(
		name 							= "overrepresentation",
		gls 							= NULL,
		glsParameterNames 				= NULL,
		transformation 					= NULL,
		transformationParameterNames 	= NULL,
		gss 							= NULL,
		gssParameterNames 				= NULL,
		globalStat 						= "global.overrepresentation",
		globalStatParameterNames 		= c("coreSet"),
		significance 					= NULL,
		significanceParameterNames 		= NULL,
		testAlternative					= NULL))
}

analysis.overrepresentation <- function(){
	return(gsAnalysis(
		name 							= "overrepresentation",
		gls 							= "gls.tStatistic",
		glsParameterNames 				= c("pValue", "alternative", "labs"),
		transformation 					= "transformation.adjustAndBinarize",
		transformationParameterNames 	= c("adjMethod", "threshold"),
		gss 							= "gss.fisherExactTest",
		gssParameterNames 				= NULL,
		globalStat 						= NULL,
		globalStatParameterNames 		= NULL,
		significance 					= NULL,
		significanceParameterNames 		= NULL,
		testAlternative					= NULL))
}
##########################################################################
#Enrichment analyses sample gls -> transformation -> gss -> significance Assessment:
#test correlation, abs.value, ar.meam and sampling
##########################################################################
analysis.averageCorrelation <- function(){
	return(gsAnalysis(
		name 							= "averageCorrelation",
		gls 							= "gls.cor",
		glsParameterNames 				= c("labs","method"),
		transformation 					= "transformation.abs",
		transformationParameterNames 	= NULL,
		gss 							= "gss.mean",
		gssParameterNames 				= NULL,
		globalStat 						= NULL,
		globalStatParameterNames 		= NULL,
		significance 					= "significance.sampling",
		significanceParameterNames 		= c("numSamples"),
		testAlternative					= "greater"))
}

##########################################################################
#Enrichment analyses sample gls -> transformation -> gss -> significance Assessment:
#test t.statistic, abs.value, ar.mean and sampling
##########################################################################
analysis.averageTStatistic <- function(){
	return(gsAnalysis(
		name 							= "averageTStatistic",
		gls 							= "gls.tStatistic",
		glsParameterNames 				= c("labs", "pValue", "alternative"),
		transformation 					= "transformation.abs",
		transformationParameterNames 	= NULL,
		gss 							= "gss.mean",
		gssParameterNames 				= NULL,
		globalStat 						= NULL,
		globalStatParameterNames 		= NULL,
		significance 					= "significance.sampling",
		significanceParameterNames 		= c("numSamples"),
		testAlternative					= "greater"))
}

##########################################################################
#global test: Goeman 2006
##########################################################################
analysis.globalTest <- function(){
	return(gsAnalysis(
		name 							= "globalTest",
		gls 							= NULL,
		glsParameterNames 				= NULL,
		transformation 					= NULL,
		transformationParameterNames 	= NULL,
		gss 							= NULL,
		gssParameterNames 				= NULL,
		globalStat 						= "global.test",
		globalStatParameterNames 		= c("response"),
		significance 					= NULL,
		significanceParameterNames 		= NULL,
		testAlternative					= NULL))
}

##########################################################################
#global ANCOVA: Hummel 2008
##########################################################################
analysis.globalAncova <- function(){
	return(gsAnalysis(
		name 							= "globalAncova",
		gls 							= NULL,
		glsParameterNames 				= NULL,
		transformation 					= NULL,
		transformationParameterNames 	= NULL,
		gss 							= NULL,
		gssParameterNames 				= NULL,
		globalStat 						= "global.ancova",
		globalStatParameterNames 		= c("group"),
		significance 					= NULL,
		significanceParameterNames 		= NULL,
		testAlternative					= NULL))
}
