#' Performs checks on the covariate.uncertainty dataframe 
#'
#' Ensures that where necessary the values are characters and that only 
#' supported sampling distributions have been selected.
#'  
#' @param covariate.uncertainty dataframe containing information used to 
#'   parametrically resample data or NULL if not required
#' @return verified covariate.uncertainty dataframe 
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @seealso \code{execute.multi.analysis}
#' @keywords input validation, data validation
#'
#'
check.covar.uncertainty <- function(covariate.uncertainty){
# check.covar.uncertainty function to perform checks on the covariate.uncertainty dataframe
#
# Arguments:  
#  covariate.uncertainty - dataframe containing information used to parametrically resample data   
#
# Value:  
#  returns the updated covariate.uncertainty dataframe
#
# Function Calls: none
# 
  if(is.null(covariate.uncertainty)){
    return(NULL)
  }            
  #Make sure these options exist and are all characters where appropriate
  character.names <- c("variable.layer", "variable.name", "cor.factor.layer", "uncertainty.layer", "uncertainty.name", "uncertainty.measure", "sampling.distribution")
  for(i in seq(along = character.names)){
    if(!character.names[i]%in%names(covariate.uncertainty)){     
      stop(paste("The ",character.names[i]," column of the covariate uncertainty dataframe has not been provided. ",sep = ""), call. = FALSE)
    }
    covariate.uncertainty[,character.names[i]] <- as.character(covariate.uncertainty[,character.names[i]])
  }
  if(!"cor.factor.name"%in%names(covariate.uncertainty)){
    stop(paste("The cor.factor.name column of the covariate uncertainty dataframe has not been provided. ",sep = ""), call. = FALSE)
  }
  if(covariate.uncertainty$cor.factor.layer == "numeric"){
    covariate.uncertainty[,"cor.factor.name"] <- as.numeric(covariate.uncertainty[,"cor.factor.name"])
  }else{
    covariate.uncertainty[,"cor.factor.name"] <- as.character(covariate.uncertainty[,"cor.factor.name"])
  }             
  #compare chosen values with those allowed
  compare <- covariate.uncertainty$sampling.distribution%in%c("Normal", "Normal.Absolute", "Lognormal.BC", "Poisson", "TruncPoisson.BC")
  if(length(which(!compare)) != 0){
    stop(paste("An unsupported sampling distribution has been chosen for covariate uncertainty. Only one of the following may be specified: Normal, Normal.Absolute, Lognormal.BC, Poisson, TruncPoisson.BC",sep = ""), call. = FALSE)
  }  
  compare <- covariate.uncertainty$variable.layer%in%c("region", "sample", "observation")
  if(length(which(!compare)) != 0){
    stop(paste("An incorrect variable layer value has been entered in the covariate uncertainty dataframe. Only one of the following may be specified: region, sample, observation.",sep = ""), call. = FALSE)
  }   
  compare <- covariate.uncertainty$uncertainty.layer%in%c("region", "sample", "observation")
  if(length(which(!compare)) != 0){
    stop(paste("An incorrect variable layer value has been entered in the covariate uncertainty dataframe. Only one of the following may be specified: region, sample, observation.",sep = ""), call. = FALSE)
  }  
  compare <- covariate.uncertainty$cor.factor.layer%in%c("numeric", "region", "sample", "observation")
  if(length(which(!compare)) != 0){
    stop(paste("An incorrect variable layer value has been entered in the covariate uncertainty dataframe. Only one of the following may be specified: numeric, region, sample, observation.",sep = ""), call. = FALSE)
  }
  compare <- covariate.uncertainty$uncertainty.measure%in%c("CV", "sd", "var")
  if(length(which(!compare)) != 0){
    stop(paste("An incorrect uncertainty measure has been entered in the covariate uncertainty dataframe. Only one of the following may be specified: CV, sd or var.",sep = ""), call. = FALSE)
  }             
  return(covariate.uncertainty)
}
