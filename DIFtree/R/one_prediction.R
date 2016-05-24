one_prediction <-
function(info,
                           item,
                           params,
                           x){
  
  node <- c() 
  all_nodes <- names(params)
  done <- FALSE
  j <- 1 
  while(!done){
    row <- info[info$number==j,]
    if(nrow(row)==1){
      var <- row$variable
      thres <- row$threshold
      if(x[var] <= thres){
        j <- node <- row$left
      } else{
        j <- node <- row$right
      }
    } else{
      done <- TRUE
    }
  }
  param_hat <- params[node==all_nodes]
  
  return(param_hat)
}
