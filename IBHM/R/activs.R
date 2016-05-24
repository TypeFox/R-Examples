# Activation functions ##############################
ActivationCandidates <- function(filter = NULL){
  # tanh activation function
  Tanh <- function(x){
    return(tanh(x));
  }
  expr(Tanh) <- 'tanh'    
  
  
  # logistic  activation function
  LogSig <- function(x){
    return(1/(1+exp(x)));
  }
  expr(LogSig) <- 'logsig'
  
  
  # linear  activation function
  Lin <- function(x){
    return(x);
  }
  expr(Lin) <- 'lin'
  
  res <- list(Tanh, LogSig, Lin)
  
  if(!is.null(filter)){
    res <- Filter(function(v) any(expr(v)==filter), res)
  }
  if(length(res)==0){
    stop('Invalid activ list given : ', paste(filter,sep=', '))
  }
  res
}
