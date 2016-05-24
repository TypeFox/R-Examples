runPC <- function(X, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                   result){
  
  # additional options for PC
  optionsList <- list("indepTest"=pcalg::gaussCItest, "fixedEdges"=NULL,
                      "NAdelete"=TRUE, "m.max"=Inf, "u2pd", 
                      "skel.method"= "stable", "conservative"=FALSE,
                      "maj.rule"=FALSE, "solve.confl"=FALSE)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  suffStat <- list(C = cor(X), n = nrow(X))
  pc.fit <- pcalg::pc(suffStat, indepTest = optionsList$indepTest, p=ncol(X), 
               alpha = alpha, 
               fixedGaps= if(is.null(variableSelMat)) NULL else (!variableSelMat), 
               fixedEdges = optionsList$fixedEdges, 
               NAdelete= optionsList$NAdelete, m.max= optionsList$m.max, 
               u2pd=optionsList$u2pd, skel.method= optionsList$skel.method, 
               conservative= optionsList$conservative, 
               maj.rule= optionsList$maj.rule, 
               solve.confl = optionsList$solve.confl, 
               verbose= verbose)
  pcmat <- as(pc.fit@graph, "matrix")
  if(directed) pcmat <- pcmat * (t(pcmat)==0)
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(as.logical(pcmat[, parentsOf[k]]))
  }
  
  result
}