
merge.setups <- function(setups){
  
  ngrp <- length(setups)
  allele.info <- NULL
  pathway <- NULL
  norm.stat <- list(score0 = list(), V = list())
  options <- setups[[1]]$options
  for(i in 1:ngrp){
    allele.info <- rbind(allele.info, setups[[i]]$allele.info)
    pathway <- rbind(pathway, setups[[i]]$pathway)
    norm.stat$score0 <- c(norm.stat$score0, setups[[i]]$norm.stat$score0)
    norm.stat$V <- c(norm.stat$V, setups[[i]]$norm.stat$V)
  }
  
  setup <- list(options = options, allele.info = allele.info, 
                pathway = pathway, norm.stat = norm.stat)
  
  setup
  
}
