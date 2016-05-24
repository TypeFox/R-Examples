fitTeachNet1 <- function(data, weights, hidden.structure, learning.rate, f, f_d,decay, m_f,er){ 
  # update weights for one hidden layer case
  # data on wich to train
  # weights for computeGrad
  # hidden.structure for number of hidden neurons
  # learning.rate in this training
  # f the acctivation function
  # f_d derivate of activation function
  # decay for update
  # m_f function to calculate  m
  # er error function
  
  H <- hidden.structure
  I <- length(weights@w_ih[,1])
  
  cost <- function(w){return((decay/2)*sum(c(sum(w@alpha^2),sum(w@alpha_h^2),sum(w@w_h^2)
                                                   ,sum(w@w_ih^2))))}
  
  # compute graddient
  grad <- Reduce("+",apply(data, 1, function(x) computeGrad1(x[2:ncol(data)],x[1],I,H,weights,f,f_d,m_f)))
  
  error_old <- er(weights,data,NA) + cost(weights)
  
  conv <- FALSE
  nichtOk <- TRUE
  t <- 1
  
  # update weights
  while(nichtOk){
    # make t so small that new weights reduce error
    t <- t*learning.rate
    weights_new <- weights - t*(grad + decay*weights)
    error <- er(weights_new,data,NA) + cost(weights_new)
    if(error<error_old){nichtOk<-FALSE}
    if(t<1e-20){nichtOk<-FALSE;print("Error can not get smaller. If you are not saticfied with the result restart procedure!");conv = TRUE}
  }

  return(list(weights_new,conv))
  
} # end of function
