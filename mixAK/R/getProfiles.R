##
##  PURPOSE:   Create a list with individual longitudinal profiles of a given variable
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   07/05/2008 (as a stand alone function)
##             26/11/2010:  added to the mixAK package
##
##  FUNCTIONS: getProfiles
##
## ==========================================================================

## *************************************************************
## getProfiles
## *************************************************************
##
getProfiles <- function(t, y, id, data)
{

## It is assumed that values for one subject follow each other
  
  T <- data[,t]
  ID <- as.character(data[,id])     ## as.character added on 20150520
                                    ## to overcome problems if id is a factor
                                    ## and data are not sorted in the same order
                                    ## as levels(id)
             ## --> this, however, created another problem...
             ## --> corrected on 20150801 (thanks to Dave Evenden)
  unID <- unique(ID)                ## added on 20150801
  ID <- factor(ID, levels = unID)   ## added on 20150801
  
  tabID <- table(ID)
  nID <- length(tabID)
  cumID <- c(0, cumsum(tabID))

  RET <- list()
  for (i in 1:nID){
    RET[[i]] <- data.frame(T[(cumID[i]+1):cumID[i+1]])
    RET[[i]] <- cbind(RET[[i]], data[(cumID[i]+1):cumID[i+1], y])
    colnames(RET[[i]]) <- c(t, y)
  }
  names(RET) <- names(tabID)
  return(RET)
}
