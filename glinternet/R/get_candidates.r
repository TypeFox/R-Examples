get_candidates = function(X, Z, res, n, pCat, pCont, numLevels, screenLimit=NULL, activeSet=NULL, norms=NULL, numCores=1){

  labels = c("cat", "cont", "catcat", "contcont", "catcont")
  candidates = list()
  candidates$variables = list()
  length(candidates$variables) = 5
  names(candidates$variables) = labels
  candidates$norms = list()
  length(candidates$norms) = 5
  names(candidates$norms) = labels
  
  #get main effect norms
  if (pCat > 0){
    candidates$variables$cat = matrix(1:pCat, ncol=1)
    if (is.null(norms$cat)) candidates$norms$cat = compute_norms_cat(X, res, n, pCat, numLevels, numCores)
    else candidates$norms$cat = norms$cat
  }
  if (pCont > 0){
    candidates$variables$cont = matrix(1:pCont, ncol=1)
    if (is.null(norms$cont)) candidates$norms$cont = compute_norms_cont(Z, res, n)
    else candidates$norms$cont = norms$cont
  }
  
  #if no limit, take all interaction pairs
  if (is.null(screenLimit)){
    if (pCat > 1) candidates$variables$catcat = cbind(rep(1:(pCat-1), (pCat-1):1), unlist(sapply(1:(pCat-1), function(x) (x+1):pCat)))
    if (pCont > 1) candidates$variables$contcont = cbind(rep(1:(pCont-1), (pCont-1):1), unlist(sapply(1:(pCont-1), function(x) (x+1):pCont)))
    if (pCat>0 && pCont>0) candidates$variables$catcont = as.matrix(expand.grid(1:pCat, 1:pCont))
  }
  else {#otherwise, sort the norms in decreasing order
    index = order(c(candidates$norms$cat, candidates$norms$cont), decreasing=TRUE)
    screenLimit = min(screenLimit, pCat+pCont)
    indices = index[1:screenLimit]
    #generate the interaction candidates
    activeCatMain = unique(c(as.vector(activeSet$catcat), activeSet$catcont[, 1]))
    catIndex = unique(c(subset(indices, indices <= pCat), activeCatMain))
    activeContMain = unique(c(as.vector(activeSet$contcont), activeSet$catcont[, 2]))
    contIndex = unique(c(subset(indices, indices > pCat)-pCat, activeContMain))
    if (pCat>1 && length(catIndex)>0) candidates$variables$catcat = cbind(rep(catIndex, (pCat-1):(pCat-length(catIndex))), unlist(sapply(1:length(catIndex), function(x) setdiff(1:pCat, catIndex[1:x]))))
    if (pCont>1 && length(contIndex)>0) candidates$variables$contcont = cbind(rep(contIndex, (pCont-1):(pCont-length(contIndex))), unlist(sapply(1:length(contIndex), function(x) setdiff(1:pCont, contIndex[1:x]))))
    if (pCat>0 && pCont>0 && length(catIndex)>0 && length(contIndex)>0) candidates$variables$catcont = as.matrix(expand.grid(catIndex, contIndex))
  }

  #get interaction norms
  if (!is.null(candidates$variables$catcat)) candidates$norms$catcat = compute_norms_cat_cat(X, res, n, numLevels, candidates$variables$catcat, numCores)
  if (!is.null(candidates$variables$contcont)) candidates$norms$contcont = compute_norms_cont_cont(Z, candidates$norms$cont, res, n, candidates$variables$contcont, numCores)
  if (!is.null(candidates$variables$catcont)) candidates$norms$catcont = compute_norms_cat_cont(X, Z, candidates$norms$cat, res, n, numLevels, candidates$variables$catcont, numCores)
  
  return(candidates)
}
