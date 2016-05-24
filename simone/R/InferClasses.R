InferStructure <- function(X, type, control) {

  clusters.crit <- control$clusters.crit
  clusters.qmin <- control$clusters.qmin
  clusters.qmax <- control$clusters.qmax
  clusters.meth <- control$clusters.meth

  ## Define setting
  if (type == "steady-state") {
    is.directed <- FALSE
  }
  if (type == "time-course") {
    is.directed <- TRUE
  }

  ## Initial Matrix Estimate
  control$penalties <- NULL # enforce NULL penalties  
  res <- simone(X, type = type, control = control)
  adj <- getNetwork(res,selection=clusters.crit)$A
  if (sum(adj) == 0) {
    cat("\nEmpty initialization network. Let me try something else...")    
    adj <- getNetwork(res)$A    
    if (sum(adj) > 0) {
      cat(" That's better, I keep this one.")    
    }
  }
  
  ## Remove Dust from analysis
  dust <- which(apply(adj,1,sum)+apply(adj,2,sum) == 0)
  ## Apply Mixer package
  if (length(dust) > 0) {
    if (length(dust) == ncol(X)) {
      stop("Empty initialization network, choose another selection criterion.")
    } else {
      classif <- mixer(adj[-dust,-dust],qmin=clusters.qmin,qmax=clusters.qmax,
                       directed=is.directed, method=clusters.meth)
    }
  } else {
    classif <- mixer(adj,qmin=clusters.qmin,qmax=clusters.qmax,
                     directed=is.directed, method=clusters.meth)
  }

  ## Retrieve best classification
  if (clusters.qmax != clusters.qmin) {
    bestCl <- classif$output[[which.max(lapply(classif$output,function(x){x$criterion}))]]
  } else {
    bestCl <- classif$output[[1]]
  }
  cl <- apply(bestCl$Taus,2,which.max)
  if (length(dust) >0) {
    clusters <- rep(0,nrow(adj))
    clusters[-dust] <- cl        
  } else {
    clusters <- cl
  }
  
  ## Weight structure
  weights <- matrix(1000,length(clusters),length(clusters))
  for (i in 1:length(clusters)) {
    for (j in 1:length(clusters)) {
      if ( !((i %in% dust)|(j %in% dust)) ) {
        weights[i,j] <- (1-bestCl$Pis[clusters[i],clusters[j]]) /max(1-bestCl$Pis)
      }
    }
  }
  
  return(list(clusters=clusters,weights=weights))
}
