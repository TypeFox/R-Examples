collapse_W <- function(W,listItems,newNames)
  {
    if(missing(newNames)) newNames <- paste("collapsedEffect",1:length(listItems),sep="")
    Wtmp1 <- W[,-unlist(listItems)]
    collapsed <- lapply(listItems, function(x) rowSums(W[,x]))
    Wtmp2 <- matrix(unlist(collapsed),ncol=length(collapsed))
    Wout <- cbind(Wtmp1,Wtmp2)
    colnames(Wout) <- c(colnames(Wtmp1),newNames)
    Wout
  }
