glmmBUGS <- function (formula, data, effects, modelFile = "model.txt", 
			initFile = "getInits.R", 
		family = c("bernoulli", "binomial", "poisson", "gaussian"), 
    spatial = NULL, spatialEffect = NULL, reparam = NULL, 
	prefix=NULL, priors=NULL,
	brugs=length(grep("unix|linux",
					.Platform$OS.type,
					ignore.case=TRUE))
)
{
    data = getDesignMatrix(formula, data, effects)
    data = na.omit(data)
    covariates = attributes(data)$covariates
    observations = attributes(data)$response
    
     #if new col created, change it back for ragged
     
     
#        reparam = paste(prefix, reparam, sep="")
     
    ragged = winBugsRaggedArray(data, effects = effects, covariates = covariates, 
        observations = observations, prefix= prefix , reparam=reparam)  
    #ragged$Xreparam = ????

	geostat=FALSE
    if (!is.null(spatial)) {
		# if spatial is a data frame, matrix, or spatialPoints or vector, use a geostatistical model
		if(class(spatial) %in% c("SpatialPoints","SpatialPointsDataFrame","complex","data.frame","matrix"))
			geostat=TRUE
	
        ragged = addSpatial(spatial, ragged, spatialEffect, prefix= prefix)
	}
	
    thepql = glmmPQLstrings(effects = effects, covariates = covariates, 
        observations = observations, data = data, family = family)
    startingValues = getStartingValues(pql = thepql, ragged = ragged, prefix=prefix, reparam = reparam )
    
    startingFunction(startingValues, file=initFile)

    spatialEffect = grep("^N[[:graph:]]+Spatial$", names(ragged), 
        value = T)
    spatialEffect = gsub("^N", "", spatialEffect)
    spatialEffect = gsub("Spatial$", "", spatialEffect)
    
    effects = paste(prefix, effects, sep="")

    
    writeBugsModel(file=modelFile, effects = effects, covariates = covariates, 
        observations = observations, family = family, spatial = spatialEffect,
		geostat=geostat,
        prefix= attributes(ragged)$prefix, reparam =reparam, priors=priors, brugs=brugs)

         return(list(ragged = ragged, startingValues = startingValues, 
        pql = thepql))
}