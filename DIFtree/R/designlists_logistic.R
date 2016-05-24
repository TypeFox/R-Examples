designlists_logistic <-
function(DM_kov,nvar,ordered_values,n_levels,n_s){
  
  v <- lapply(1:nvar,function(j) 1:(n_levels[j]-1)) 
  w <- lapply(1:nvar, function(j) rep(paste0("s",j),n_s[j]))
  
  design_upper <- lapply(1:nvar, function(j){
    design_matrix <- sapply(1:(n_levels[j]-1),function(k) { ifelse(DM_kov[,j] > ordered_values[[j]][k],1,0)})
    colnames(design_matrix) <- paste0(w[[j]],v[[j]],"_u")
    design_matrix
  })
  design_lower <- lapply(1:nvar, function(j){
    design_matrix <- abs(design_upper[[j]]-1)
    colnames(design_matrix) <- paste0(w[[j]],v[[j]],"_l")
    design_matrix
  })
  return(list(design_lower,design_upper))
}
