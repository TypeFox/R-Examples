tree_Logistic <-
function(Y,
                          DM_kov,
                          npersons,
                          nitems,
                          nvar,
                          sumscore,
                          ordered_values,
                          n_levels,
                          n_s,
                          type,
                          alpha,
                          nperm,
                          trace){
  
  # call functions 
  if(type=="udif"){
    to_return <- tree_udif(Y,DM_kov,npersons,nitems,nvar,sumscore,ordered_values,n_levels,n_s,alpha,nperm,trace)
  }
  if(type=="dif"){
    to_return <- tree_dif(Y,DM_kov,npersons,nitems,nvar,sumscore,ordered_values,n_levels,n_s,alpha,nperm,trace)
  }
  if(type=="nudif"){
    to_return <- tree_nudif(Y,DM_kov,npersons,nitems,nvar,sumscore,ordered_values,n_levels,n_s,alpha,nperm,trace)
  }
  
  return(to_return)

}
