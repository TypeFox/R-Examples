runGES <- function(X, parentsOf, variableSelMat, setOptions, directed, verbose, 
                   result){
  
  # additional options for GES
  optionsList <- list("turning"=TRUE, "maxDegree"=integer(0))
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  vers <- unlist(packageVersion('pcalg'))[1:2]
  score <- new("GaussL0penObsScore", X)
  vers <- unlist(packageVersion('pcalg'))[1:2]
  tryNewVersion <- try(
    {
      G <- pcalg::ges(score, 
               fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
               turning = optionsList$turning, maxDegree=optionsList$maxDegree, 
               verbose=verbose)
    }, 
    silent=TRUE)
  if(class(tryNewVersion)=="try-error"){
    G <- pcalg::ges(ncol(X), score, 
             fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
             turning = optionsList$turning, maxDegree=optionsList$maxDegree, 
             verbose=verbose)
  }
  gesmat <- as(G$essgraph, "matrix")
  if(directed) gesmat <- gesmat * (t(gesmat)==0)
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(gesmat[, parentsOf[k]] == 1) 
  }
  
  result
}