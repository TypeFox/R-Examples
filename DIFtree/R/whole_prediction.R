whole_prediction <-
function(info,
                             item,
                             params,
                             X){
  
  params_hat <- c()
  n_pred     <- nrow(X)
  for(i in 1:n_pred){
    params_hat[i] <- one_prediction(info,item,params,X[i,])
  }
  return(params_hat)
}
