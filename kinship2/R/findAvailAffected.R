# Automatically generated from all.nw using noweb

findAvailAffected <- function(ped, avail, affstatus)
  ## Try trimming one subject by affection status indicator
  ## If ties for bits removed, randomly select one of the subjects

  {
 
    notParent <- !is.parent(ped$id, ped$findex, ped$mindex)
    
    if(is.na(affstatus)) {
      possiblyTrim <- ped$id[notParent & avail & is.na(ped$affected)]
    } else {
      possiblyTrim <- ped$id[notParent & avail & ped$affected==affstatus]
    }
    nTrim <- length(possiblyTrim)
    
    if(nTrim == 0)
      {
        return(list(ped=ped,
                    idTrimmed = NA,
                    isTrimmed = FALSE,
                    bitSize = bitSize(ped)$bitSize))
      }
    
    trimDat <- NULL
    
  for(idTrim in possiblyTrim) {

    
      avail.try <- avail
      avail.try[ped$id==idTrim] <- FALSE
      id.rm <- findUnavailable(ped, avail.try)
      newPed <- pedigree.trim(id.rm, ped)
      trimDat <- rbind(trimDat,
                  c(id=idTrim, bitSize=bitSize(newPed)$bitSize))
    }

    bits <- trimDat[,2]

    # trim by subject with min bits. This trims fewer subject than
    # using max(bits).

    idTrim <- trimDat[bits==min(bits), 1]
    
    ## break ties by random choice
    if(length(idTrim) > 1)
      {
        rord <- order(runif(length(idTrim)))
        idTrim <- idTrim[rord][1]
      }

    
    avail[ped$id==idTrim] <- FALSE
    id.rm <- findUnavailable(ped, avail)
    newPed <- pedigree.trim(id.rm, ped)
    pedSize <- bitSize(newPed)$bitSize
    avail <- avail[!(ped$id %in% id.rm)]

    return(list(ped=newPed,
                newAvail = avail,
                idTrimmed = idTrim,
                isTrimmed = TRUE,
                bitSize = pedSize))
  }

