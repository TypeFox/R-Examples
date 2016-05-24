as.party.quint <-
function(quint.out, nodeID=1L, ...){
  quintNodes <- to.party(nodeID,quint.out,...)
  party.object <- party(quintNodes,quint.out$data
                        , fitted=data.frame(
                            "(fitted)"=colnames(quint.out$nind)[
                              apply(quint.out$nind==1,1,which)],
                            "(response)"=quint.out$data[,1],
                            check.names=FALSE)
                        , terms=terms(as.formula(paste(
                            colnames(quint.out$data)[1],"~",
                            paste(colnames(quint.out$data)[-1],collapse="+"))
                            )),
                        ...
                        )
  party.object$ni <- quint.out$li
  class(party.object) <- c("constparty",class(party.object))
  return(party.object)
}
