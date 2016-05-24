getBetaWeights <-
function(model, documents=NULL) {
  logbeta <- model$beta$logbeta
  K <- model$settings$dim$K
  vocab <- model$vocab
  
  #Let's start by marginalizing
  margbeta <- exp(logbeta[[1]])
  if(length(logbeta) > 1) {
    weights <- model$settings$covariates$betaindex
    tab <- table(weights)
    weights <- tab/sum(tab)
    #marginalize
    margbeta <- margbeta*weights[1]
    for(i in 2:length(model$beta$logbeta)) {
      margbeta <- margbeta + exp(model$beta$logbeta[[i]])*weights[i]
    }
  }
  
  ##
  # figure out how to weight the topics.
  # NB: if they didn't provide topics use naive weights
  #     otherwise calibrate thetas by the total counts
  #     per document.
  if(is.null(documents)) {
    weights <- colSums(model$theta)
  } else {
    D.n <- unlist(lapply(documents, function(x) sum(x[2,])))
    weights <- colSums(D.n*model$theta)
  }
  
  return(list(beta=margbeta, weights=weights))
}
