runRFCI <- function(X, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                    result){
  
  # additional options for RFCI
  optionsList <- list("indepTest"=pcalg::gaussCItest,
                      "skel.method"="stable", "fixedEdges"=NULL,
                      "NAdelete"=TRUE, "m.max"=Inf,"rules"=rep(TRUE,10),
                      "conservative"=FALSE, "maj.rule"=FALSE)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  suffStat <- list(C = cor(X), n = nrow(X))
  rfci.fit <- pcalg::rfci(suffStat, indepTest = optionsList$indepTest, 
                   p=ncol(X), alpha=alpha, 
                   fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
                   fixedEdges=optionsList$fixedEdges, 
                   NAdelete=optionsList$NAdelete, m.max=optionsList$m.max, 
                   skel.method= optionsList$skel.method, 
                   conservative= optionsList$conservative, 
                   maj.rule=optionsList$maj.rule, rules=optionsList$rules, 
                   verbose= verbose )
  rfcimat <- as(rfci.fit@amat, "matrix")
  if(directed) rfcimat <- rfcimat * (t(rfcimat)==0)
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(as.logical(rfcimat[, parentsOf[k]]))
  }

  result
}