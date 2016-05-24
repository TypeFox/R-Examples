####################################################################################
####################################################################################
## FUNCTION TO CONVERT MICROSIMULATION OUTPUT INTO WIDE FORMAT                    ##
## SZ, April 2014                                                                 ##
####################################################################################
####################################################################################

convertToWideFormat <- function(pop){
  
  giveSeq <- function(nu){
    return(1:nu)
  }
  
  popTemp <- pop
  ns <-  data.frame(table(popTemp$ID),stringsAsFactors=FALSE)
  ns <- ns[order(as.numeric(as.character(ns[,1]))),]
  colnames(ns) <- c('ID','ns')
  popTemp <- merge(popTemp, ns, by='ID')
  popTemp <- popTemp[order(as.numeric(popTemp[,c('ID')])),] 
  nsU <- popTemp$ns[which(!duplicated(popTemp[,'ID']))]
  popTemp$Episode <- unlist(sapply(nsU,giveSeq))
  popTemp <- popTemp[,c('ID','birthDate','initState','ns','Episode','From','To','transitionTime','transitionAge')]
  popWide <- reshape(popTemp, timevar = 'Episode', idvar = 'ID', direction = 'wide', 
                     v.names=c('From','To', 'transitionTime', 'transitionAge'))  
  popWide$ns[which(is.na(popWide$transitionTime.1))] <- 0
  return(popWide)
}

