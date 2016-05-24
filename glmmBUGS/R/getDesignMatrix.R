getCovList = function(covariates, data, effects) {

  if(is.list(covariates)) covariates = unique(unlist(covariates))

  # find out which covariates apply at each level
  covList = list()

  # get effects factors which are unique 
  # each ID only appears once, regardless of the parent level  
  uniqueEffects=as.data.frame(matrix(NA, nrow=dim(data)[1],
    ncol=length(effects), dimnames=list(NULL, effects)))

 
  theeffect = uniqueEffects[,1] = data[,effects[1]]
  for(Deffect in effects[-1]) {
    uniqueEffects[[Deffect]] = theeffect = paste(theeffect, data[,Deffect], sep="/")
  }

  # loop through covariates to see which level they apply at
  Ndata = dim(data)[1]
# return(list(a=uniqueEffects, b=data$P0))
  for(Dcov in covariates) {
    Deffect = length(effects)
    keepgoing=TRUE
    while( keepgoing & (Deffect>0)) {
      # if there is only one value of this covariate for each observation with 
      #  this effect
      
      theorder = order(uniqueEffects[,effects[Deffect]], data[,Dcov])
      # put the data in order
      data2 = data[theorder,Dcov]
      # and the effects in the same order
      effect2 = uniqueEffects[theorder,effects[Deffect]]
      # find the places where the data changes values
      datadiff = which(data2[-1] != data2[-Ndata])
      # and where the effect changes
      effectdiff = which(effect2[-1] != effect2[-Ndata])

      # if the effect always changes when the data changes, the data must be 
      #  constant within an effect
      keepgoing = all(datadiff %in% effectdiff)

    
      Deffect = Deffect - keepgoing  
    }
    Deffect = Deffect + 1
    if(Deffect <= length(effects))
      covList[[effects[Deffect] ]] = c( covList[[effects[Deffect] ]], Dcov)
    else
      covList[["observations"]] = c(covList[["observations"]], Dcov)  
  }
  
  return(covList)
}
 
 
 getDesignMatrix = function(formula, data, effects=NULL) {

  covariates = attr(terms(formula), "term.labels")
  interactions = grep(":", covariates, value=TRUE)
  mainEffects = covariates[! covariates %in% interactions]
  
  response = as.character(attributes(terms(formula))$variables)[2]
  response = unlist(strsplit(response, "\\+"))
  response=gsub("[[:space:]]", "", response)

  if(!all(mainEffects %in% names(data)))
    warning(mainEffects[!mainEffects %in% names(data)], " not found")

  # find which covariates apply at which level
  covList = getCovList(mainEffects, data, effects)
  
  # make the most populous category the baseline
  bases = getBases(covList, data, effects)

  designMatrix = model.matrix(formula, data, contrast.arg=bases)
  # get rid of intercept column
  designMatrix = designMatrix[,dimnames(designMatrix)[[2]] != "(Intercept)", drop=FALSE]

  covariates = dimnames(designMatrix)[[2]]
  # replace : with _
  covariates = gsub(":", "_", covariates)
  dimnames(designMatrix)[[2]] = covariates
  
  stuff = merge(data[,response, drop=FALSE], designMatrix,
    by="row.names", all=T)
  data = merge(stuff, data[,effects, drop=FALSE], by.x="Row.names",
    by.y="row.names", all=T)
  rownames(data) = data[,"Row.names"]

  covList = getCovList(covariates, data, effects)
  
  attr(data, "covariates") = covList
  attr(data, "response") = response

	if(is.logical(data[[response[1]]]))
		data[[response]] = as.numeric(data[[response[1]]])
	if(is.factor(data[[response[1]]]))
		warning("response can't be specified as a factor")
	if(is.character(data[[response[1]]]))
		warning("response can't be specified as a character string")


  return(data)  
  
}

getBases = function(covariates, data,  effects=NULL) {
    if(!is.list(covariates))
      covariates = list(observations = covariates)

    # find the reference category of factors
    bases = list()  
    
    allfactors = names(which(sapply(data, is.factor)))

 
    for(D in covariates$observations[covariates$observations %in% allfactors])  {
      # remove spaces from labels
      levels(data[[D]]) <- gsub("[[:space:]]", "", levels(data[[D]]))
    
      bases[[D]] = contr.treatment(levels(data[[D]]), 
    order(table(data[[D]]), decreasing=TRUE)[1])
    
   } 


    if(is.null(effects)) {
      effects = names(covariates)[names(covariates) != "observations"]
    }


    if(length(effects)) {
    data2 = data
    for(Deffect in length(effects):1) {

       data2 = data2[!duplicated(data2[,effects[seq(Deffect, 1)] ]),]

       thecov =  covariates[[effects[Deffect]]]

       thecov = thecov[thecov %in% allfactors]

       for(D in thecov)
             bases[[D]] = contr.treatment(levels(data[[D]]), 
                order(table(data[[D]]), decreasing=TRUE)[1])
    }
    }
     bases

}
