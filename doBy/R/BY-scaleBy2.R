scaleBy <- function(formula, data, center=TRUE, scale=TRUE, details=0){
  vv <- .get_variables(formula, data=data, id=NULL, debug.info=1)
  lhs.num <- vv$lhs.num
  rhs.grp <- vv$rhs.grp
  if (details>0){
    cat(sprintf("grouping variables = %s\n", toString(rhs.grp)))
    cat(sprintf("scale variables    = %s\n", toString(lhs.num)))
  }
  #str(vv)
  rh.string <- .get_rhs_string(data, rhs.grp)
  rh.unique <- unique(rh.string)
  #print(rh.string)
  for (uu in rh.unique){
    idx <- uu == rh.string
    #print(which(idx))
    data[idx, lhs.num ] <- scale(data[idx, lhs.num, drop=FALSE], 
                                 center=center, scale=scale)
  }	  
  data
}
