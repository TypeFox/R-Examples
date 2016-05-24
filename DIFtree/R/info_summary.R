info_summary <-
function(splits,
                         item,
                         model,
                         type){
  
  if(model==2 & type==2){
    dif_items <- unique(c(splits[[1]][,"item"],splits[[2]][,"item"]))
  } else{
    dif_items <- unique(splits[,"item"])
  }
  
  dif <- ifelse(item %in% dif_items,"yes","no")
  if(dif=="yes"){
    if(model==1 | (model==2 & type==1)){
      type      <- "uniform"
      variables <- paste(unique(splits[splits[,"item"]==item,"variable"]),collapse=",")
      nos    <- nrow(splits[splits[,"item"]==item,])
    } 
    if(model==2 & type==2){
      type      <- ifelse(item %in% unique(splits[[2]][,"item"]),"non-uniform","uniform")
      variables <- paste(unique(c(splits[[1]][splits[[1]][,"item"]==item,"variable"],splits[[2]][splits[[2]][,"item"]==item,"variable"])),collapse=",")
      nos       <- nrow(rbind(splits[[1]][splits[[1]][,"item"]==item,],splits[[2]][splits[[2]][,"item"]==item,]))
    }
    if(model==2 & type==3){
      type      <- "non-uniform"
      variables <- paste(unique(splits[splits[,"item"]==item,"variable"]),collapse=",")
      nos       <- nrow(splits[splits[,"item"]==item,])
    }
  } else{
    type <- variables <- nos <- "---"
    
  }

  
  output <- list("item"=item,
                 "dif"=dif,
                 "type"=type,
                 "variables"=variables,
                 "nos"=nos)
  return(output)
}
