##########################################################################
#enrichment analysis: general procedure
##########################################################################
doGLS <- function(
		dat,
		analysis,
		parameters){

	if(!is.null(analysis$gls)){
		glsParameters <- getNeededParameters(formals(analysis$gls),
			parameters[analysis$glsParameterNames])
		gls <- do.call(analysis$gls, c(list(dat), glsParameters))
	}else{
		gls <- dat
	}
	return(gls)
}

doTransformation <- function(
		analysis,
		parameters,
		gls){

	if(!is.null(analysis$transformation)){
		transformationParameters <- getNeededParameters(formals(analysis$transformation),
			parameters[analysis$transformationParameterNames])
		transformation <- do.call(analysis$transformation, c(list(gls), transformationParameters))
	}else{
		transformation <- gls
	}
	return(transformation)
}

doGSS <- function(
		dat,
		geneSet,
		analysis,
		parameters,
		transformation){

	geneSetIndices <- which(tolower(rownames(dat)) %in% tolower(geneSet))
	if(!is.null(analysis$gss)){
		gssParameters <- getNeededParameters(formals(analysis$gss),
			parameters[analysis$gssParameterNames])
		gss <- do.call(analysis$gss, c(list(transformation, geneSetIndices = geneSetIndices), gssParameters))
	}else{
		gss <- transformation
	}
	return(gss)
}

doGlobalAnalysis <- function(
		dat,
		geneSet,
		analysis,
		parameters){

	globalParameters <- getNeededParameters(formals(analysis$globalStat),
		parameters[analysis$globalStatParameterNames])
	gs <- do.call(analysis$globalStat,
		c(parameters[!(names(parameters)%in%names(globalParameters))],
			list(dat = dat, geneSet = geneSet), globalParameters))
	return(gs)
}

doSignificanceAssessment <- function(
		dat,
		geneSet,
		analysis,
		glsValues,
		parameters){

	if(is.null(analysis$significance)){
		gsSig <- NULL
	}else{
		sigParameters <- getNeededParameters(formals(analysis$significance),
			parameters[analysis$significanceParameterNames])

		gsSig <- do.call(analysis$significance,
			c(parameters[!(names(parameters)%in%names(sigParameters))],
				list(dat = dat,
				geneSet = geneSet,
				analysis = analysis,
				glsValues = glsValues), sigParameters))
	}
	return(gsSig)
}
