prep.data.event.AJ <-
function(dNs, Ys, sum_dNs,states, tr.states) {
  
  split_dNs <- strsplit(colnames(dNs), "  ") ## string splits names
  split_Ys <- strsplit(colnames(Ys), " ")  ## string split names of Ys
  split_sum_dNs <- strsplit(colnames(sum_dNs), " ") ## string splits names of dNs
  
  ## reducing dNs & Ys to just event times & noncens states
  ##  looks at noncensored columns
  event.dNs.id <- which(sapply(split_dNs, function(x) x[1]=="tr 1 2" | x[1]=="tr 1 3" | x[1]=="tr 2 3"))
  ## identifies times where transitions occur
  event.row.id <- which(apply(dNs[, event.dNs.id, drop=FALSE], 1, function(x) any(x>0)))
  dNs.event <- dNs[event.row.id, event.dNs.id, drop=FALSE] ## reduces dNs  
  
  tr.states.event <- names(which(sapply(tr.states, function(x) length(x)>0))) 
  event.Ys.id <- which(sapply(split_Ys, function(x) x[2]%in%tr.states.event)) 
  Ys.event <- Ys[event.row.id, event.Ys.id, drop=FALSE] ## reduces Ys
  
  event.sum_dNs.id <- which(sapply(split_sum_dNs, function(x) x[2]%in%states))
  sum_dNs.event <- sum_dNs[event.row.id, event.sum_dNs.id, drop=FALSE]
  
  ans <- list(dNs=dNs.event, Ys=Ys.event, sum_dNs=sum_dNs.event)
  return(ans)
  
}
