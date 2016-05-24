fitTeachNet2 <- function(data, weights, hidden.structure, learning.rate, f, f_d,decay, m_f,er){ 
  # update weights for one hidden layer case
  
  # data on wich to train
  # weights for computeGrad
  # hidden.structure for number of hidden neurons
  # learning.rate for this training
  # f the acctivation function
  # f_d derivate of activation function
  # decay for update
  # m_f function to calculate m
  # er error function
  
  H <- hidden.structure[2]
  M <- hidden.structure[1]
  I <- length(weights@w_im[,1])
  
  cost <- function(w){return((decay/2)*sum(c(sum(w@alpha^2),sum(w@alpha_2h^2),sum(w@alpha_1m^2),
                                                   sum(w@w_h^2),sum(w@q_mh^2), sum(w@w_im^2))))}
  
  # compute graddient
  grad <- Reduce("+",apply(data, 1, function(x) computeGrad2(x[2:ncol(data)],x[1],I,M,H,weights,f,f_d,m_f)))
  
  error_old <- er(weights,data,H) + cost(weights)
  
  conv <- FALSE
  nichtOk <- TRUE 
  t <- 1
  
  # update weights
  while(nichtOk){
    t <- t*learning.rate
    weights_new <- weights - t*(grad + decay*weights)
    error <- er(weights_new,data,H) + + cost(weights_new)
    if(error<error_old){nichtOk<-FALSE}
    if(t<1e-20){nichtOk<-FALSE;print("Error can not get smaller. If you are not saticfied with the result restart procedure!");conv = TRUE}
  }
  
  return(list(weights_new,conv))
  
} # end of function