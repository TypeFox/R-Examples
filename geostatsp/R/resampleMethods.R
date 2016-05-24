setGeneric('resampleMethods', function(
				formula, covariates, data=NULL) 
			standardGeneric("resampleMethods"))

setMethod("resampleMethods", 
    signature("ANY", "ANY", "missing"), 
    function(formula, covariates, data=NULL){
      data = data.frame()
      callGeneric(formula, covariates, data)
    }
)


setMethod("resampleMethods", 
    signature("ANY", "ANY", "NULL"), 
    function(formula, covariates, data=NULL){
      data = data.frame()
      callGeneric(formula, covariates, data)
    }
)

setMethod("resampleMethods", 
    signature("ANY", "ANY", "SpatialPointsDataFrame"), 
    function(formula, covariates, data=NULL){
      data=data@data
      callGeneric(formula, covariates, data)
    }
)


# convert covariates to a list      
setMethod("resampleMethods", 
    signature("ANY", "Raster", "data.frame"), 
    function(formula, covariates, data=NULL){
      covariatesList = vector('list', nlayers(covariates))
      names(covariatesList) = names(covariates)
      for(D in names(covariates))
        covariatesList[[D]] = covariates[[D]]
      covariates = covariatesList
      callGeneric(formula, covariates, data)
    }
)


setMethod("resampleMethods", 
    signature("character", "list", "ANY"), 
    function(formula, covariates, data=NULL){

      # restrict covariates to those listed in formula
      covariates = covariates[
          intersect(formula, names(covariates))]
      # ignore formula
      formula = ~1
      callGeneric(formula, covariates, data)
    }
)



setMethod("resampleMethods", 
		signature("formula", "list", "data.frame"), 
		function(formula, covariates, data=NULL){
# decide which method to use when reprojecting covariates
			# factors must be ngb, numerics are bilinear
      
      allVars = all.vars(formula)
      allVars = intersect(allVars, names(covariates))
      
      allterms =rownames(attributes(terms(
                  update.formula(formula, junk~.)))$factors)
      factorsInFormula = grep("^factor\\(", allterms, value=TRUE)
      factorsInFormula = gsub("^factor\\(|\\)$", "", factorsInFormula)
      
      factorsInCovariates = unlist(lapply(covariates, is.factor))
      factorsInCovariates=names(factorsInCovariates)[factorsInCovariates]

      varsInData = intersect(allVars, names(data))
      factorsInData = unlist(
          lapply(data[,varsInData, drop=FALSE], is.factor)
      )
      factorsInData = names(factorsInData)[factorsInData]

      method = rep("bilinear", length(names(covariates)))
      names(method)=names(covariates)
      method[names(method) %in% 
              c(factorsInFormula, factorsInCovariates, factorsInData)
      ] = "ngb" 
			method
    }
)

			