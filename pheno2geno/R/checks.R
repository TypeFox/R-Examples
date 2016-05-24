#
# checks.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Apr, 2013
# first written Dec, 2011
# Contains: numericCheck.internal, genotypeCheck.internal, inRangeCheck.internal 
#           inListCheck.internal, check.population, crossContainsMap.internal 
#           defaultCheck.internal
#

# genotypeCheck.internal
#
# DESCRIPTION:
#  checking if given object is containing only 0,1 and NAs
# OUTPUT:
#  boolean
#
genotypeCheck.internal <- function(objectToBeChecked, genotypes, allow.na=FALSE){
  converted <- as.numeric(as.matrix(objectToBeChecked))
  if(any(is.na(converted))){
    if(!(allow.na)) return(FALSE)
    if(sum(is.na(converted))==sum(is.na(objectToBeChecked))){
      nrOfCorrect <- sum(is.na(converted))
      for(genotype in genotypes){ nrOfCorrect <- nrOfCorrect + sum(converted==genotype,na.rm=TRUE) }
      if(nrOfCorrect==length(converted)) return(TRUE)
      return(FALSE)
    }
    return(FALSE)
  }
  return(TRUE)
}

# numericCheck.internal
#
# DESCRIPTION:
#  Checking if given object is numeric or could be converted to numeric
# OUTPUT:
#  boolean
#
numericCheck.internal <- function(objectToBeChecked, allow.na=FALSE){
  converted <- as.numeric(as.matrix(objectToBeChecked))
  if(any(is.na(converted))){
    if(!(allow.na)) return(FALSE)
      if(sum(is.na(converted)) == sum(is.na((objectToBeChecked)))) return(TRUE)
    return(FALSE)
  }
  return(TRUE)
}

############################################################################################################
#                                          *** check.population ***
#
# DESCRIPTION:
#  checking if given object is a correct population object
# OUTPUT:
#  none
############################################################################################################
check.population <- function(x,verbose=FALSE){
  #if(length(class(x))!=2) stop("Incorrect class of the object.\n")
  if(class(x)[1]!="population") stop("Object is not of a class population.\n")
  if(!(class(x)[2] %in% c("riself", "f2", "bc", "risib"))) stop("Type of the population: ",class(x)[2]," not recognized.\n")
  if(is.null(x$founders$phenotypes) && !("noParents" %in% x$flags) && !("annots" %in% x$flags)) stop("No founders phenotype data found, this is not a valid object of class population.\n")
  if(is.null(x$offspring$phenotypes)) stop("No offspring phenotype data found, this is not a valid object of class population.\n")
  if(is.null(x$founders$groups) && !("noParents" %in% x$flags) && !("annots" %in% x$flags)) stop("No founders groups information found, this is not a valid object of class population.\n")
}

############################################################################################################
#                  *** defaultCheck.internal  ***
#
# DESCRIPTION:
#   making sure that default parameter is used, when parameter is speicified by =c("","")
# OUTPUT:
#  default parameter from list of possible
############################################################################################################
defaultCheck.internal <- function(parameterToBeChecked,nameOfParameter,maxLength,defVal){
  if(length(parameterToBeChecked) == maxLength){
    invisible(defVal)
  }else if(length(parameterToBeChecked) != 1){
    stop("wrong parameter ",nameOfParameter," length, choose one out of possible\n")
  }else{
    invisible(parameterToBeChecked)
  }
}

