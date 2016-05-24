export.nodal <-
function(blauObj, niches = TRUE){
  if (is.null(blauObj$isInNiche)){
    print("Nothing to export.")
  }

  if (niches == TRUE){
    if ("ecologyNames" %in% colnames(blauObj$isInNiche)){
      to.export <- cbind(blauObj$ids,blauObj$isInNiche[, 1:(ncol(blauObj$isInNiche)-1)])
    }
    else{
      to.export <- cbind(blauObj$ids,blauObj$isInNiche)
    }
  }
  else{
    to.export <- data.frame(matrix(0, nrow = nrow(blauObj$nodalLocal), ncol=0))
  }

  if (!is.null(blauObj$nodalLocal)){
    to.export <- cbind(to.export, blauObj$nodalLocal)
  }
  if (!is.null(blauObj$nodalGlobal)){
    to.export <- cbind(to.export, blauObj$nodalGlobal)
  }
  if (!is.null(blauObj$nodalNetwork)){
    to.export <- cbind(to.export, blauObj$nodalNetwork)
  }
  rownames(to.export) <- NULL
  return(to.export)
}
