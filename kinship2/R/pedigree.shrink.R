# Automatically generated from all.nw using noweb

#$Log: pedigree.shrink.q,v $
#Revision 1.5  2010/09/03 21:11:16  sinnwell
#add shrunk "avail" vector to result, keep status and affected in pedObj
#
#Revision 1.4  2010/09/03 19:15:03  sinnwell
#add avail arg which is not part of ped object.  Re-make ped object at the end with status and affected, if given
#
#Revision 1.2  2009/11/17 23:08:18  sinnwell
#*** empty log message ***
#
#Revision 1.1  2008/07/16 20:23:07  sinnwell
#Initial revision
#
pedigree.shrink <- function(ped, avail, affected=NULL, seed=NULL, maxBits = 16){
  if(class(ped) != "pedigree")
    stop("Must be a pegigree object.\n")
   
  ## set the seed for random selections
  if(is.null(seed))
    {
      seed <- sample(2^20, size=1)
    }
  set.seed(seed)

  if(any(is.na(avail)))
    stop("NA values not allowed in avail vector.")
  
  if(is.null(affected))
    affected = if(is.matrix(ped$affected)) ped$affected[,1] else ped$affected

  ped$affected = affected
 
 
  idTrimmed <- numeric()
  idList <- list()
  nOriginal <- length(ped$id)
 
  bitSizeOriginal <- bitSize(ped)$bitSize
  
  ## first find unavailable subjects to remove anyone who is not 
  ## available and does not have an available descendant
  
  idTrimUnavail <- findUnavailable(ped, avail)

  
  if(length(idTrimUnavail)) {    
    
    pedTrimmed <- pedigree.trim(idTrimUnavail, ped)
    avail <- avail[match(pedTrimmed$id, ped$id)]
    idTrimmed <- c(idTrimmed, idTrimUnavail)
    idList$unavail <- paste(idTrimUnavail, collapse=' ')

  } else {
    ## no trimming, reset to original ped
    pedTrimmed <- ped
  }

  
  ## Next trim any available terminal subjects with unknown phenotype
  ## but only if both parents are available
  
  ## added nNew>0 check because no need to trim anymore if empty ped
  
  nChange <- 1
  idList$noninform = NULL
  nNew <- length(pedTrimmed$id)

  while(nChange > 0 & nNew > 0){
    nOld <- length(pedTrimmed$id)
    
    ## findAvailNonInform finds non-informative, but after suggesting 
    ## their removal, checks for more unavailable subjects before returning
    idTrimNonInform <- findAvailNonInform(pedTrimmed, avail)
    
    if(length(idTrimNonInform)) {
        pedNew <- pedigree.trim(idTrimNonInform, pedTrimmed)
        avail <- avail[match(pedNew$id, pedTrimmed$id)]
        idTrimmed <- c(idTrimmed, idTrimNonInform)
        idList$noninform = paste(c(idList$noninform, 
               idTrimNonInform), collapse=' ')
        pedTrimmed <- pedNew
        
    }
    nNew <- length(pedTrimmed$id)
    nChange <- nOld - nNew
    
  }
  
  ##  Determine number of subjects & bitSize after initial trimming
  nIntermed <- length(pedTrimmed$id)
  
  bitSize <- bitSize(pedTrimmed)$bitSize
    
  ## Now sequentially shrink to fit bitSize <= maxBits
    
  bitVec <- c(bitSizeOriginal,bitSize)
  
  isTrimmed <- TRUE
  idList$affect=NULL 
  
  while(isTrimmed & (bitSize > maxBits))
    {  
        
      ## First, try trimming by unknown status
      save <- findAvailAffected(pedTrimmed, avail, affstatus=NA)
      isTrimmed <- save$isTrimmed
      
      ## Second, try trimming by unaffected status if no unknowns to trim
      if(!isTrimmed)
        {
          save <- findAvailAffected(pedTrimmed, avail, affstatus=0)
          isTrimmed <- save$isTrimmed
          
        }
      
      
      ## Third, try trimming by affected status if no unknowns & no unaffecteds
      ## to trim
      if(!isTrimmed) {
        save <- findAvailAffected(pedTrimmed, avail, affstatus=1)
        isTrimmed <- save$isTrimmed
      }
      
      if(isTrimmed)  {
        pedTrimmed <- save$ped
        avail <- save$newAvail
        bitSize <- save$bitSize
        bitVec <- c(bitVec, bitSize)          
        idTrimmed <- c(idTrimmed, save$idTrimmed)
        idList$affect = paste(c(idList$affect, save$idTrimmed), 
          collapse=' ')
      }
      
      
    } # end while (isTrimmed) & (bitSize > maxBits)
  
  
  nFinal <- length(pedTrimmed$id)
  
  obj <- list(pedObj = pedTrimmed,
              idTrimmed = idTrimmed,
              idList = idList,
              bitSize = bitVec,
              avail=avail,
              pedSizeOriginal = nOriginal,
              pedSizeIntermed = nIntermed,
              pedSizeFinal  = nFinal,
              seed = seed)


  oldClass(obj) <- "pedigree.shrink"

  return(obj)
} 


