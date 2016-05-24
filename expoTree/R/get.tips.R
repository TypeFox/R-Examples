get.tips <- function(node.ages,code) {
  off <- which(node.ages[,1]==code)
  off.codes <- node.ages[off,2]
  tips <- lapply(off.codes,get.tips,node.ages=node.ages)
  return(rbind(cbind(off,off.codes),do.call(rbind,tips)))
}

