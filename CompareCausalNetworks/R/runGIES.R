runGIES <- function(X, interventions, parentsOf, variableSelMat, setOptions, 
                    directed, verbose, result){

    # check validity of input arguments for GIES
  if(is.null(interventions)) 
    stop("'interventions' cannot be 'NULL' for method 'gies'")
  
  # additional optionsList for GIES
  optionsList <- list("turning"=TRUE, "maxDegree"=integer(0))
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)

  isn <- which(sapply(interventions, is.null))
  if(length(isn)>0) { for (k in isn) interventions[[k]] <- numeric(0) }
  targets <- unique(interventions)
  target.index <- match(interventions, targets)
  targets <- as.list(lapply(targets, as.integer))
  
  if(any(table(unlist(targets)) == length(targets)))
    stop(paste("One variable is intervened on in all settings. At least one\n",
         "setting needs to be entirely different for 'gies' to run."), 
         call. = FALSE)
  
  score <- new("GaussL0penIntScore", data=X, targets=targets,  
               target.index = target.index)
  
  tryNewVersion <- try(
    {
      tmp <- pcalg::gies(
        score, 
        fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
        turning=optionsList$turning, maxDegree=optionsList$maxDegree,
        verbose=verbose)
    },
    silent=TRUE)
  
  if(class(tryNewVersion)=="try-error"){
    tmp <- pcalg::gies(
      ncol(X), targets, score, 
      fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat),
      turning=optionsList$turning, maxDegree=optionsList$maxDegree,
      verbose=verbose)
  }
  giesmat <- as(tmp$essgraph, "matrix")
  if(directed) giesmat <- giesmat * (t(giesmat)==0)
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(giesmat[, parentsOf[k]] == 1) 
  }
  
  result
}