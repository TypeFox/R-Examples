print.nodeHarvest <-
function(x,nonodes=3,...){


  Z <- x[["nodes"]]
  
  cat("\n", "\t Node Harvest estimator:", "\n \t", length(Z), "selected nodes"  )
  cat("\n","\t Nodes contain, in the weighted average,", signif(length(x[["Y"]])/sum(sapply(Z,attr,"weight")),3),"observations")
  cat("\n \t Node means range from", signif(min(sapply(Z,attr,"mean")),3),"to",signif(max(sapply(Z,attr,"mean")),3))
  cat("\n","\t Estimator explains", signif(100*(1- mean( (x[["Y"]]-x[["predicted"]])^2)/mean( (x[["Y"]]-mean(x[["Y"]]))^2) ),3),"% of the variance on training data")

  nonodes <- min(nonodes,length(Z))
  cat("\n \t", "The", nonodes, "nodes with largest weight:")
  weights <- sapply(Z,attr,"weight")
  ord <- order(weights,decreasing=TRUE)
  tmp <- 0
  for (kk in ord[1:nonodes]){
    tmp <- tmp+1
    cat("\n\n \t ",tmp,") Node ",kk,", containing ",attr(Z[[kk]],"n")," training observations, with mean ", signif(attr(Z[[kk]],"mean"),3)," and weight ", signif( attr(Z[[kk]],"weight"),3)," :",sep="")
    if(attr(Z[[kk]],"depth")>0){
      node <- drawtext(Z,kk,varnames=x[["varnames"]],plot=FALSE)
      for (lc in 1:length(node)){
        cat("\n\t\t ", node[[lc]])
      }
    }else{
      cat("\n\t\t ROOT NODE")
    }
  }
  cat("\n")
    
  
  
  
   
}

