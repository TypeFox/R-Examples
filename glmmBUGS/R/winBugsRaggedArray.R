`winBugsRaggedArray` <-
function(data, 
    effects = names(data)[-length(names(data))], 
    covariates=NULL, 
    observations=names(data)[length(names(data))],
    returnData=FALSE, prefix=NULL, reparam=FALSE) {

  # if a vector of covariates is specified, assume they're all observation level
  #if(!is.list(covariates)) covariates = list(observations = covariates)
  if(!is.list(covariates)) covariates = list(covariates)
  
  # check to see covariate categories correspond to effects
  covNames= ! (names(covariates) %in% effects)
  if(any(covNames)) {
     if(sum(covNames>1)){#add here
      warning("names of covariates don't match effect names")
     names(covariates)[covNames] = "observations" 
     }#added here
  }
   
 
  result <- list()
  Neffects <- length(effects)

  # convert everything to characters to make sure it's all in the correct order
  for(D in effects) {
    dataD = as.character(data[,D])
    data[,D] = encodeString(dataD,max(nchar(dataD)))
  }

  # reorder the data
  if(length(effects)>1) {theorder = do.call(order, data[,effects])
  }else theorder = order(data[,effects])
  data=data[theorder,]

  result[[paste("N", prefix, effects[1], sep="")]] = length(unique(data[,effects[1]]))

  # get the ragged array sequences
  if(Neffects>1) {
    for(D in 2:Neffects) {
      data[,effects[D]] = paste(data[,effects[D-1]], data[,effects[D]], sep="/")  
      result[[paste("S", prefix, effects[D-1], sep="")]] = 
        getRaggedSeq(data[,effects[c(D-1, D)]])
    }
  }
  Sfull = seq(1, dim(data)[1])

  # observations
  result[[paste("S", prefix, effects[Neffects], sep="")]] =
    getRaggedSeq(cbind(data[,effects[Neffects]], Sfull)) 
  for(Dobservation in observations)
    result[[paste(prefix, Dobservation, sep="")]] = data[,Dobservation] 

   # get covariates
 

  #Dlevel = "observations"
   Dlevel = paste(prefix,effects[1],sep="")

      
  # change data dataframe to a matrix
  data = as.matrix(data[,unlist(covariates), drop=FALSE])

  if(!is.null(covariates[[Dlevel]])) {
      result[[paste("X",  Dlevel, sep="")]] =
        data[Sfull,covariates[[Dlevel]]] 
  }
  
    # obseration level covariates
  Plevel = paste(prefix, "observations", sep="")
    
  if(!is.null(covariates[[Plevel]])) {
      result[[paste("X",  Plevel, sep="")]] =
        data[Sfull,covariates[[Plevel]]] 
  }
  

  # the other levels  
  for(Dlevel in rev(effects)) {
     Dlevel = paste(prefix, Dlevel, sep="")
    # extract the covariates at this level
    Sfull = Sfull[result[[paste("S", Dlevel, sep="")]]]
    Sfull = Sfull[-length(Sfull)] 
    if(!is.null(covariates[[Dlevel]])) {
      theXname =  paste("X", Dlevel, sep="")
      theSname =  paste("S", Dlevel, sep="")
      result[[theXname]] =  
        data[Sfull,covariates[[Dlevel]]] 
    # give row names to the covariates
      if(is.matrix(result[[theXname]])) { 
        rownames( result[[theXname]]) =
          names(result[[theSname]])[-length(result[[theSname]])]
      } else {
        names( result[[theXname]]) =
          names(result[[theSname]])[-length(result[[theSname]])]
      }   
    }
  }
  

 # add reparam 
 
 if(is.character(reparam))   {
reparamName = reparam
reparam = list()
for(D in reparamName)
reparam[[D]] = list()
 } 
 
for(D in names(reparam)){
     if(D %in% names(covariates)){
       theXname= paste("X", D, "reparam", sep="")  
       if(length(reparam[[D]])){
         result[[theXname]] = reparam[[D]]
         }else{result[[theXname]] = apply(as.matrix(data[, covariates[[D]]]), 2, mean) }
       
      } 
 }    
     

attributes(result)$prefix = prefix
attributes(result)$covariates = covariates

 if(returnData) 
    result=list(data=data, result=result)
  
  
    return(result)
}

