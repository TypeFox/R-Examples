sae.train <- function(x,hidden=c(10),
                      activationfun="sigm",
                      learningrate=0.8,
                      momentum=0.5,
                      learningrate_scale=1,
                      output="sigm",
                      numepochs=3,batchsize=100,
                      hidden_dropout=0,visible_dropout=0.2                      
                      ){
  if (!is.matrix(x)) 
    stop("x must be a matrix!")
  input_dim <- ncol(x)
  size <- c(input_dim, hidden)
  sae <- list(
    input_dim = input_dim,
    hidden = hidden,
    size = size
  )
  train_x <- x
  message(sprintf("training layer 1 autoencoder ..."))
  sae$encoder[[1]] <-  nn.train(train_x,train_x,hidden=c(hidden[1]),
                                      activationfun=activationfun,
                                      learningrate=learningrate,
                                      momentum=momentum,
                                      learningrate_scale=learningrate_scale,
                                      output=output,
                                      numepochs=numepochs,batchsize=batchsize,
                                      hidden_dropout=hidden_dropout,visible_dropout=visible_dropout)
  
  if(length(sae$size) > 2){
    for(i in 2:(length(sae$size) - 1)){
      pre <- t( sae$encoder[[i-1]]$W[[1]] %*% t(train_x) + sae$encoder[[i-1]]$B[[i-1]] )
      if(sae$encoder[[i-1]]$activationfun == "sigm"){
        post <- sigm( pre )
      }else if(sae$encoder[[i-1]]$activationfun == "tanh"){
        post <- tanh(pre)
      }else{
        stop("unsupport activation function 'nn$activationfun'");
      }  
      train_x <- post
      message(sprintf("training layer %d autoencoder ...",i))
      sae$encoder[[i]] <- nn.train(train_x,train_x,hidden=c(hidden[i]),
                                   activationfun=activationfun,
                                   learningrate=learningrate,
                                   momentum=momentum,
                                   learningrate_scale=learningrate_scale,
                                   output=output,
                                   numepochs=numepochs,batchsize=batchsize,
                                   hidden_dropout=hidden_dropout,visible_dropout=visible_dropout)
    }
  }
  sae
}


