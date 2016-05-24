#' Checks the list of model names supplied by the user
#'
#' Performs various checks on the models and model names supplied by the user, 
#' including checking that the data in each model are the same across species 
#' codes to ensure valid comparison using the AIC/AICc/BIC selection criteria.
#'
#' In addition, this routine establishes whether the analyses are double
#' observer or not and whether the observations are individuals of clusters.
#'
#' @param model.names a list of vectors of model names. Each list element is
#'   named according to the species code. 
#' @param ddf.models a list of ddf objects named in model.names
#' @return list of two boolean values relating to cluster size and double
#'   observer analyses. 
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @seealso \code{execute.multi.analysis}
#' @keywords input validation, data validation
#'

check.ddf.models <- function(model.names, ddf.models){             
# 
# check.ddf.models - Performs various checks on the model names supplied by the user
#
# Arguments:
#
#  model.names   - a list of vectors specifying the model names for each species code
#  ddf.models    - a list of ddf objects named in model.names   
#
# Value: list of two boolean values relating to cluster size and double observer analyses. 
# 
# Function calls: 
#   is.same - local function to compare two datasets to see if they match
#                                      

#add a check that all the detected are 1 in the ddf data [used to be a problem in mrds]
  
  is.same <- function(data1, data2){
  # is.same function to compare two datasets to see if they match
  #
  # Arguments:
  #   data1 - dataframe for comparison
  #   data2 - dataframe for comparison
  #
  # Value: returns numeric
  #   0 - dataframes are the same
  #   1 - differ in number of observations
  #   2 - differ in object ID and/or distances
  #  
    #if the dataframes have different lengths return 1
    if(nrow(data1) != nrow(data2)){
      return(1)
    }
    #If they are of the same dimension check if they are identical
    if(ncol(data1) == ncol(data2)){
      compare <- data1 == data2
      #if the dataframes are identical return 0
      if(length(which(compare == FALSE)) == 0){
        return(0)
      }
    }
    #order data by object and compare object id's and distances
    compare.data1 <- data1[order(data1$object), c("object","distance")]
    compare.data2 <- data2[order(data2$object), c("object","distance")]
    compare <- compare.data1 == compare.data2
    if(length(which(compare == FALSE)) == 0){
      return(0)
    }else {
      return(2)
    }
  }
  

  species.name <- names(model.names)
  clusters <- FALSE
  double.observer <- NULL
  
  #CHECK WHETHER ALL MODELS ARE PROVIDED AND IF IT IS A MR DOUBLE OBSERVER ANALYSIS 
  model.type <- NULL
  counter <- 1
  for(sp in seq(along = model.names)){     
    #for every model
    for(m in seq(along = model.names[[sp]])){
      #method <- try(ddf.models[[model.names[[sp]][m]]]$method, silent = TRUE)
      method <- ddf.models[[model.names[[sp]][m]]]$method
      #CHECK MODEL EXISTS
      if(is.null(method)){
        #ddf object doesn't exist    
        stop(paste("ddf object ",m,", analysis name ",model.names[[sp]][m],", for species code ",species.name[sp]," has not been provided.",sep = ""), call. = FALSE)
      } 
      model.type[counter] <- method 
      counter <- counter + 1
    }#next model
  }# next species
  double.observer <- which(model.type%in%c("trial", "trial.fi", "io", "io.fi"))
  ds <- which(model.type%in%c("ds"))
  #unsupported <- which(!model.type%in%c("trial", "trial.fi", "io", "io.fi", "ds"))
  unsupported <- which(!model.type%in%c("ds"))
  if(length(unsupported) > 0){
    stop(paste("Unsupported model types have been selected: ",paste(model.type[unsupported], collapse = ", "), sep = ""), call. = FALSE)
  }
  if(length(double.observer) == length(model.type)){
    double.observer <- TRUE
    #check all are trial or all are io
    if(!(length(which(model.type%in%c("trial", "trial.fi"))) == length(model.type) | length(which(model.type%in%c("io", "io.fi"))) == length(model.type))){  
      stop(paste("Models must either be all trial or all io, not a mixture.", sep = ""), call. = FALSE)
    }
  }else if(length(ds) == length(model.type)){
    double.observer <- FALSE
  }else if(length(double.observer) > 0 & length(ds) > 0){
    stop("Models must either be all double observer mark-recapture or all standard distance sampling models, not a mixture.", call. = FALSE)
  }
  rm(model.type, counter, ds, unsupported)

  #CHECK THAT MODELS FOR EACH SPECIES ARE UNIQUE
  for(sp in seq(along = model.names)){
    assoc.model.names <- model.names[[sp]]
    for(m in seq(along = assoc.model.names)){
      for(mcheck in seq(along = assoc.model.names)){
        if(m == mcheck){
          next
        }else if(assoc.model.names[m] == assoc.model.names[mcheck]){      
          stop(paste("The model names are not unique for species ",names(model.names)[sp],".", sep = ""), call. = FALSE)
        }
      }#next model for checking
    }#next model   
  }#next species
  #CHECK DATA MATCHES ACROSS DIFFERENT MODELS FOR THE SAME SPECIES
  for(sp in seq(along = model.names)){     
    for(m in seq(along = model.names[[sp]])){ 
      #check model exists 
      ddf.data <- ddf.models[[model.names[[sp]][m]]]$data
      if(m == 1){
        #Get first dataset to compare all others too
        check.data <- ddf.data
      }else if(is.same(check.data, ddf.data) != 0){
        stop(paste("Datasets within species must contain the same data to ensure the model selection criteria are valid. The ",species.name[sp]," analyses ",model.names[[sp]][1]," and ",model.names[[sp]][m]," do not have the same sightings and/or associated distances", sep = ""), call. = FALSE)
      }
    }#next model
    #CHECK IF DATA CONTAINS CLUSTER SIZES "size" (either all must or all most not)
    if("size"%in%names(ddf.data)){
      if(sp > 1 & !clusters){
        stop("Cluster size must be present in all datasets within the ddf models or none.", call. = FALSE)
      }
      clusters <- TRUE
    }else if(clusters){
      stop("Cluster size must be present in all datasets within the ddf models or none.", call. = FALSE)
    }    
  }#next species
  return(list(clusters = clusters, double.observer = double.observer))
}