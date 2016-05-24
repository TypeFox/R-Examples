createBlockIDs <- function(obj, data, id.var){
  
  if(length(obj[[1]]) == 1){
    obj.simp <- obj[[1]]$`1`
  }else{
    obj.simp <- NULL
    for(gp in obj[[1]]){
      obj.simp <- rbind(obj.simp, gp)
    }
  }
  
  row.n <- nrow(obj.simp)
  bbb <- rep(NA, nrow(data)) 
  
  for(col.idx in 1:(ncol(obj.simp)-1)){ ## only for level.two == FALSE
    tmp.colname <- paste("col", col.idx, sep = "")
    assign(tmp.colname, obj.simp[, col.idx])   
    tmp.col <- get(tmp.colname)
    
    if(is.factor(tmp.col)){
      tmp.col <- unfactor(tmp.col)
    }
    
    for(row.idx in 1:length(tmp.col)){
      bbb[data[[id.var]] == tmp.col[row.idx]] <- row.idx
    }
  }
  
  return(bbb)
}
