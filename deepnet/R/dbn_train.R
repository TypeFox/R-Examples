dbn.train <- function(x,hidden=c(10,10),
                      numepochs=3,batchsize=100,
                      learningrate=0.8,learningrate_scale=1,momentum=0.5,
                      visible_type="bin",hidden_type="bin",cd=1){
  if (!is.matrix(x)) 
    stop("x must be a matrix!")
  input_dim <- ncol(x)  
  dbn <- list(
    size = c(input_dim, hidden)
  )
  train_x <- x
  message("training layer 1 rbm ...")
  dbn$rbm[[1]] <- rbm.train(train_x,hidden[1],
                            numepochs=numepochs,batchsize=batchsize,
                            learningrate=learningrate,learningrate_scale=learningrate_scale,
                            momentum=momentum,
                            visible_type=visible_type,hidden_type=hidden_type,cd=cd)
  
  if(length(dbn$size) > 2){
    for(i in 2:(length(dbn$size) - 1)){
      train_x <- rbm.up(dbn$rbm[[i-1]], train_x)  
      message(sprintf("training layer %d rbm ...",i))
      dbn$rbm[[i]] <- rbm.train(train_x,hidden[i],
                                numepochs=numepochs,batchsize=batchsize,
                                learningrate=learningrate,learningrate_scale=learningrate_scale,
                                momentum=momentum,
                                visible_type=visible_type,hidden_type=hidden_type,cd=cd)
    }
  }
  dbn
}

dbn.down <- function(dbn,h,round=10){
  hi <- h
  i <- length(dbn$size) - 1 #top rbm
  for(j in 1:round){
    vi <- rbm.down(dbn$rbm[[i]],hi)
    hi <- rbm.up(dbn$rbm[[i]],vi)
  }
  if(length(dbn$size) > 2){
    hi <- vi
    for(i in (length(dbn$size) - 2):1){
      vi <- rbm.down(dbn$rbm[[i]],hi)
      hi <- vi
    }    
  }
  vi
}
