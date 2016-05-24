`getItemDetails` <- function(raschResult, item, class=1, camelCase=TRUE){
  if(substr(raschResult$model,1,3)=="mix"){
    rR <- raschResult[[1]][[class]]$item.par
  } else{
    rR <- raschResult$item.par
  }
  
  itemNames <- colnames(raschResult$item.par$delta)
  if(item  %in% itemNames){
    itemSelect <- itemNames == item	
  } else itemSelect <- item
  
  n.cat <- rR$tau[,itemSelect]
  n.cat <- sum(! is.na(n.cat))
  if(is.null(rR$SE.tau)){
    noTau <- TRUE 
  } else if(is.na(rR$SE.tau)){
    noTau <- TRUE
  } else{
    noTau <- FALSE
  }
  if(noTau){
    SE.tau <- NA
  } else{
    SE.tau <- rR$SE.tau[1:n.cat,itemSelect]
  }	
  out <- list(item.name=item, n.cat=n.cat+1, delta.i = rR$delta.i[itemSelect], SE.delta.i = rR$SE.delta.i[itemSelect], 
              tau = rR$tau[1:n.cat,itemSelect], SE.tau = SE.tau, infit = rR$in.out[,"infit"][itemSelect], in.Z = rR$in.out[,"in.Z"][itemSelect],
			  outfit = rR$in.out[,"outfit"][itemSelect],out.Z = rR$in.out[,"out.Z"][itemSelect],
			  itemMean = rR$itemDescriptives[itemSelect,"itemMean"], pBis = rR$itemDescriptives[itemSelect,"pBis"], bis = rR$itemDescriptives[itemSelect,"bis"])	
  
  if(camelCase) names(out) <- c("itemName", "nCat", "deltaI", "SEDeltaI", "tau", "SETau", "infit", 
                                "inZ", "outfit", "outZ", "itemMean", "pBis", "bis")
  out
}