# 
# Copyright (c) 2010, 2014, IBM Corp. All rights reserved. 
# 		
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>. 
# 
# 

idaListModels <- function() {
  idaQuery("call idax.list_models()")
}

idaDropModel <- function(modelname) {
  xx <- parseTableName(modelname);
  model <- paste('"',xx$schema,'"."',xx$table,'"',sep=''); 	
  
  idaQuery("CALL IDAX.DROP_MODEL('model=",model,"')")
}

idaGetModelname <- function(object) {
  if(inherits(object,c("rules", "sequencerules"))) {
    return(object@info$model)
  } else if(inherits(object,c("idaKMeans","idaNaiveBayes","idaLm"))){
    return(object$model)
  } else {
    stop("Model type not supported")    
  }
}

idaRetrieveModel <- function(modelname) {
  xx <- parseTableName(modelname);
  
  model <- xx$table
  modelSchema <- xx$schema
  
  models <- idaListModels()
  modelAlgorithm <- models$ALGORITHM[model == models$MODELNAME &
                                     modelSchema == models$MODELSCHEMA]
  if(length(modelAlgorithm)==0){
    stop("The model you are trying to retrieve does not exist.")
  }
  #THE NAME OF THE RETRIEVEMODEL FUNCTION HAS TO BE AT THE SAME INDEX AS THE NAME OF THE
  #RESPECTIVE ALGORITHM...
  retrievemethods <- c("idaRetrieveKMeansModel", "idaRetrieveRulesModel", "idaRetrieveNBModel",
                     "idaRetrieveSeqRulesModel", "idaRetrieveidaLmModel", "idaRetrieveTreeModel",
                     "idaRetrieveTreeModel")
  
  algorithms <- c("Kmeans", "Association Rules", "Naive Bayes", "Sequential Patterns",
                "Linear Regression", "Regression Tree", "Decision Tree")
  
 
    algNr <- match(modelAlgorithm, algorithms)
    if(is.na(algNr)){stop(paste("The algorithm", modelAlgorithm, "is not supported by ibmdbR."))}
    tryCatch(
    model <- get(retrievemethods[algNr])(modelname),
    error = function(e){
             stop(e)
           }
  )
  return(model)
}

idaModelExists <- function(modelname) {
  xx <- parseTableName(modelname);
  models <- idaListModels()
  return(nrow(models[(models$MODELNAME==xx$table)&(models$MODELSCHEMA==xx$schema),])>0)
}

idaGetValidModelName <- function(prefix) {
  prefix <-toupper(prefix); 
  while (TRUE) {
    name <- paste(prefix, floor(runif(1,0,100000)),'_',as.integer(Sys.time()), sep="")
    if (!idaModelExists(name))
      return(name)
  }
}