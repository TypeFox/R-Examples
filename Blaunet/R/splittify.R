splittify <-
function(blauObj, ecologyId, ecologyRows) {
  miniBlau <- list() #this object holds the ROW parts of the blau object that are in the relevant ecology
  miniBlau$ids <- blauObj$ids[ecologyRows, , drop=FALSE]
  miniBlau$memberships <- blauObj$memberships[ecologyRows, , drop=FALSE]
  miniBlau$weights <- blauObj$weights[ecologyRows, , drop=FALSE]
  miniBlau$dimensions <- blauObj$dimensions[ecologyRows, , drop=FALSE]

  if (!is.null(blauObj$primMemCol)){
    miniBlau$primMemCol <- blauObj$primMemCol
  }
  if (!is.null(blauObj$isInNiche)){
    miniBlau$isInNiche <- blauObj$isInNiche[which(blauObj$isInNiche[,'ecologyNames'] == ecologyId), which(colnames(blauObj$isInNiche) != 'ecologyNames'), drop=FALSE] #this picks the rows in the ecology, and all columns except for the one with the ecology names
  }

  if (!is.null(blauObj$topbounds) && !is.null(blauObj$lowbounds)){
    miniBlau$topbounds <- blauObj$topbounds[which(blauObj$topbounds[,'ecologyNames'] == ecologyId), which(colnames(blauObj$topbounds) != 'ecologyNames'), drop=FALSE] 
    miniBlau$lowbounds <- blauObj$lowbounds[which(blauObj$lowbounds[,'ecologyNames'] == ecologyId), which(colnames(blauObj$lowbounds) != 'ecologyNames'), drop=FALSE]     
  }

  #keep the full network because the networked function will delete extraneous vertices
  miniBlau$graph <- blauObj$graph

  return(miniBlau)
}
