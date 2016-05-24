ape2time <- function(tree=NULL,eps=0.0) {
  node.ages <- node.age(tree)
  # branching events
  bev <- node.ages[,4]>0 & node.ages[,5]>1
  # sampling events
  sev <- node.ages[,4]==0 & node.ages[,5]>0
  btimes <- as.numeric(node.ages[bev,3])
  stimes <- as.numeric(node.ages[sev,3])
  times <- c(btimes,stimes)
  ttypes <- as.integer(c(btimes*0+1,stimes*0))
  times <- max(times)-times
  extant <- times <= eps
  times <- times[!extant]
  ttypes <- ttypes[!extant]
  root <- which(is.na(node.ages[,1]))
  if (length(root) > 0) {
    times <- c(times,max(times)+node.ages[root,3]-node.ages[root,6])
    ttypes <- c(ttypes,1)
  }
  o <- order(times)
  cbind(times[o],ttypes[o])
}

