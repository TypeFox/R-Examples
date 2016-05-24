
plot.marginalRelevance <-  function (x, newdata = NULL,  n.feat = NULL,  type = "default", ...)
{
  
  NF <- dim(x$orderedData)[2]
  
  if(is.null(n.feat))n.feat <- NF 
  
  if(is.null(newdata)){
    NTR <- dim(x$orderedData)[1]
    ordered <- x$orderedData
  }
  else{
    NTR <- dim(newdata)[1]
    ordered <- matrix(0, NTR, NF) 
    for (i in 1 : NF)ordered[ ,i] <- newdata[ ,which(x$rank == i)]
  }
  
  #if(is.null(col))col <- 1 # rep(1, NTR)
  if (type == "parallelcoord"){
    
    dt <- ordered[, 1 : n.feat]
    rdt <- apply(dt, 2L, range, na.rm = TRUE)
    dt <- apply(dt, 2L, function(dt) (dt - min(dt, na.rm = TRUE))/(max(dt, na.rm = TRUE) - min(dt, na.rm = TRUE)))
    
    matplot(1L:ncol(dt), t(dt), xaxt = "n", yaxt = "n", bty="n",  xlab = "Ordered by marginal relevance",  ylab ="", type = 'l',lty = 1, ...) 
    
    axis(1, at = 1 : n.feat, labels = x$bestVars[1 : n.feat], cex.axis = 0.7)
    for (i in 1L:ncol(dt)) lines(c(i, i), c(0, 1), col = "grey70")
    
  } 
  else if (type == "pairs"){
    pairs(ordered[, 1 : n.feat], labels = x$bestVars[1 : n.feat], ...)
  }  
  else if (type == "default"){
    plot(1 : NF, x$score, xlab = "",  ylab = "MR score" ,  col = 0)
    lines(1 : NF, x$score)
  }
  else  stop("error:  Plot type not supported for a marginalRelevance object")
}
