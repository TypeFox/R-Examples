to.party <-
function(nodeID,quint.out, ...){
  if(class(nodeID)=="quint"){quint.out <- nodeID; nodeID <- 1L}
  if(is.null(quint.out$var.names)) quint.out$var.names <- colnames(quint.out$data)
  if(nodeID %in% quint.out$li[,1]) return(partynode(id=as.integer(nodeID),
                                                   info = quint.out$li[quint.out$li[,1]==nodeID,]
                                                   ))
  if(nodeID %in% quint.out$si[,1]) return(partynode(id=as.integer(nodeID),
                                                   split =  partysplit(varid = as.integer(
                                                                         which(
                                                                               quint.out$var.names==quint.out$si[quint.out$si[,1]==nodeID,3])),
                                                     breaks = quint.out$si[quint.out$si[,1]==nodeID,4]),
                                                   ##index if factor
                                                   kids = lapply(c(nodeID*2,nodeID*2+1),to.party,quint.out),
                                                   )
                                         )
}
