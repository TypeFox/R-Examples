`summary.msc` <-
function(object, ...){
  num.Crits <- length(object$msc.List)
  num.Sumfuns <- length(object$summary.Functions)
  corners = list()
  for(i in 1:num.Crits){corners[[i]]<-object$Sum.Stats[which.max(object$Sum.Stats[[i]]),]}
  names(corners) <- names(object$msc.List)
  minima = list()
  for(j in 1:num.Sumfuns){minima[[j]]<-object$Sum.Stats[order(object$Sum.Stats[,num.Crits+j])[1:5],]}
  names(minima) = names(object$summary.Functions)
  temp <- list(Corners=corners, Minima=minima)
  class(temp) <- "summary.msc"
  temp
}

