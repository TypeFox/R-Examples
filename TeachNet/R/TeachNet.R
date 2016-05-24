TeachNet <- function(data, hidden.structure=-1, threshold.TotalError=0.1, stepMax=100,
                     learning.rate=0.9, acc.fct="logistic", err.fct="sse",startWeights=NULL, decay=0, 
                     sameSample=FALSE,sampleLength=0.7,all=FALSE,eval=TRUE){ 
  # trains the neural network
  
  is.data(data)                        # check for data
  is.numberOfNeurons(hidden.structure) # check for correct hidden.structure 
  is.thres.error(threshold.TotalError) # check for correct threshold.TotalError
  is.stepMax(stepMax)                  # check for correct stepMax
  is.learn(learning.rate)              # check for correct learning.rate
  is.acct(acc.fct)                     # check for correct acc.fct
  is.err(err.fct)                      # check for correct err.fct
  is.decay(decay)                      # check for correct decay
  is.sample(sameSample)                # check for correct sameSample
  is.sampleLeng(sampleLength)          # check for correct err.fct
  is.sample(all)                       # check for correct all
  is.sample(eval)                      # check for correct eval

  
  data[2:ncol(data)] <- scale(data[2:ncol(data)])
  
  if(!all){
    # build train and test data 
    if(!sameSample){
      num <- as.integer(nrow(data)*sampleLength)
      train <- data[sample(1:nrow(data), num), ]
      test <- data
    } else {
      good <- subset(data, data[1]==0)
      bad <- subset(data, data[1]==1)
      length.good <- nrow(good)
      length.bad <- nrow(bad)
      if(length.bad<length.good){    
        length.sample <- round(sampleLength*length.bad)
        ende <- length.bad
      } else {length.sample <- round(sampleLength*length.good);ende <- length.good}
      
      good <- good[sample(length.good), ]
      bad <- bad[sample(length.bad), ]
      bad.train <- bad[1:length.sample, ]
      bad.test <- bad[(length.sample+1):length.bad, ]
      good.train <- good[1:length.sample, ]
      good.test <- good[(length.sample+1):ende, ]
      
      train <- rbind(good.train, bad.train)
      train <- train[sample(nrow(train)), ]
      
      test <- rbind(good.test, bad.test)
      test <- test[sample(nrow(test)), ]
    }
  } else {
    train <- data
    test <- data
  }
  
  # count variables
  numberOfVar <- ncol(data)
  
  # initialize some values
  steps <- 0
  isOk <- TRUE
  
  # set acctivation function
  if(acc.fct=="logistic"){ f <- function(x){return(logistic(x))}; f_d <- function(x){return(logistic.differential(x))}  }
  
  # set error function
  if(err.fct=="sse"){ e <- function(weights,data,h2){return(sumSquaredError(weights,data,h2))}; 
                      m_f <-function(z,y){return(2*( f(z) - y ))} }
  if(err.fct=="ce"){ e <- function(weights,data,h2){return(sumCrossEntropy(weights,data,h2))}; 
                     m_f <-function(z,y){return((-1)*((y/f(z))-((1-y)/(1-f(z)))))} }
  
  # test for automatic rule
  if(hidden.structure[1]==-1) {hidden.structure=(as.integer(numberOfVar/2))}
  
  # set function for hidden layer
  if(is.na(hidden.structure[2])) {
    fitTeachNet <- function(train,weights, hidden.structure=hidden.structure, learning.rate=learning.rate,f=f,
                          f_d=f_d,decay=decay,m_f=m_f,er=e){return(fitTeachNet1(train,weights, hidden.structure=hidden.structure, 
                                                                              learning.rate=learning.rate,f=f, f_d=f_d,decay=decay,m_f=m_f,er=er))}
    
    createWeights <- function(I,hidden.structure){return(createWeights1(ncol(data)-1,hidden.structure[1]))}
    
    computeOutput <- function(x,weights){return( computeOutput1(x,weights) )}
    
    costWeights <- function(){return((decay/2)*sum(c(sum(weights@alpha^2),sum(weights@alpha_h^2),sum(weights@w_h^2)
                                                     ,sum(weights@w_ih^2))))}
    
  } else {
    
    fitTeachNet <- function(train,weights, hidden.structure=hidden.structure, learning.rate=learning.rate,f=f,
                          f_d=f_d,decay=decay,m_f=m_f,er=e){return(fitTeachNet2(train,weights, hidden.structure=hidden.structure, 
                                                                              learning.rate=learning.rate,f=f, f_d=f_d,decay=decay,m_f=m_f,er=er))}
    
    createWeights <- function(I,hidden.structure){return(createWeights2((ncol(data)-1),hidden.structure))}
    
    computeOutput <- function(x,weights){return( computeOutput2(x,weights) )}
    
    costWeights <- function(){return((decay/2)*sum(c(sum(weights@alpha^2),sum(weights@alpha_2h^2),sum(weights@alpha_1m^2),
                                                     sum(weights@w_h^2),sum(weights@q_mh^2), sum(weights@w_im^2))))}
    
  }
  
  # set weights
  if (is.null(startWeights)) {
    weights <- createWeights(ncol(data)-1,hidden.structure)
  } else {weights=startWeights}
  
  # while loop to repeat dataset
  while(isOk) {
    
    # calculate next weights   
    
    ret <- fitTeachNet(train,weights, hidden.structure=hidden.structure, learning.rate=learning.rate,f=f,
                     f_d=f_d,decay=decay,m_f=m_f)
    
    weights <- ret[[1]]
    
    steps <- steps + 1
    
    # print some information every 10 steps and check for reached total error
    if (steps%%1==0){
      cat("Number of processed Updates = ", steps , ",")
      error <- e(weights,train,hidden.structure[2])
      error <- error + costWeights()
      cat("train error = ", error , "\n")
      # end if train error is small enough
      if(error < threshold.TotalError){isOk=FALSE}
    }
    
    if(ret[[2]]){cat("Reached steps = ", steps , "\n");isOk=FALSE}
    
    # end if stepmax is reached
    if(steps>=stepMax){isOk=FALSE}
    
  }# end of while
  
  if(steps < 10){error <- e(weights,train,hidden.structure[2])+ costWeights()}
  
  cat("Reached train error = ", error , "\n")
  
  # evaluation
  if(eval){
    # predict for test data
    pred <- predict(weights, test)
      
    # evaluate on test data
    first <- find.Threshold(test[,1], stepsize=0.01, predict=pred)
    thresh <- first$Threshold[1]
    conf <- confusion(pred,test[,1], threshold=thresh)
    acc <- accuracy.me(test[,1],predict=pred,thres=thresh)
      
    print("Evaluation:")
    print(first)
    print("Confusion Matrix:")
    print(conf)
    print("Modell accuracy:")
    print(acc)
  }
  
  return(weights)
  
}# end of function