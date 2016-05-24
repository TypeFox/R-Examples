##' Training Neural Network
##'
##' Training single or mutiple hidden layers neural network by BP 
##' @param x matrix of x values for examples
##' @param y vector or matrix of target values for examples
##' @param initW initial weights. If missing chosen at random 
##' @param initB initial bias. If missing chosen at random 
##' @param hidden vector for number of units of hidden layers.Default is c(10).
##' @param activationfun activation function of hidden unit.Can be "sigm","linear" or "tanh".Default is "sigm" for logistic function
##' @param learningrate learning rate for gradient descent. Default is 0.8.
##' @param momentum momentum for gradient descent. Default is 0.5 .
##' @param learningrate_scale  learning rate will be mutiplied by this scale after every iteration. Default is 1 .
##' @param numepochs number of iteration for samples  Default is 3.
##' @param batchsize size of mini-batch. Default is 100.
##' @param output function of output unit, can be "sigm","linear" or "softmax". Default is "sigm".
##' @param hidden_dropout drop out fraction for hidden layer. Default is 0.
##' @param visible_dropout drop out fraction for input layer Default is 0.
##' @examples
##' Var1 <- c(rnorm(50,1,0.5),rnorm(50,-0.6,0.2))
##' Var2 <- c(rnorm(50,-0.8,0.2),rnorm(50,2,1))
##' x <- matrix(c(Var1,Var2),nrow=100,ncol=2)
##' y <- c(rep(1,50),rep(0,50))
##' nn <-nn.train(x,y,hidden=c(5))
##' @author Xiao Rong
##' @export
nn.train <- function(x,y,initW=NULL,initB=NULL,hidden=c(10),
                    activationfun="sigm",
                    learningrate=0.8,
                    momentum=0.5,
                    learningrate_scale=1,
                    output="sigm",
										numepochs=3,batchsize=100,
                    hidden_dropout=0,visible_dropout=0) {
	if (!is.matrix(x)) 
     stop("x must be a matrix!")
	input_dim <- ncol(x)
	output_dim <- 0;
	if(is.vector(y)){
		output_dim <- 1
	}else if(is.matrix(y)){
		output_dim <- ncol(y)
	}
	if (output_dim == 0) 
        stop("y must be a vector or matrix!")
	size <- c(input_dim, hidden, output_dim)
	vW <- list() 
  vB <- list()
  if(is.null(initW) || is.null(initB)){
	  W <- list()
  	B <- list()
	  #random init weights and bias between layers							 
	  for( i in 2:length(size) ){
		  W[[i-1]] <- matrix(runif(size[i]*size[i-1],min=-0.1,max=0.1), c(size[i],size[i-1]));
		  B[[i-1]] <- runif(size[i],min=--0.1,max=0.1);
		  vW[[i-1]] <- matrix(rep(0,size[i]*size[i-1]),c(size[i],size[i-1]))
      vB[[i-1]] <- rep(0,size[i])
	  }
  }else{
    W <- initW
    B <- initB
    for( i in 2:length(size) ){
      vW[[i-1]] <- matrix(rep(0,size[i]*size[i-1]),c(size[i],size[i-1]))
      vB[[i-1]] <- rep(0,size[i])
      if(nrow(W[[i-1]]) != size[i] || ncol(W[[i-1]]) != size[i-1] ){
        stop("init W size is not eq to network size!")  
      }    
      if(length(B[[i-1]]) != size[i]){
        stop("init B size is not eq to network size!")  
      }
    }
  }
	
	nn <- list(
		input_dim = input_dim,
		output_dim = output_dim,
		hidden = hidden,
		size = size,
		activationfun = activationfun,
		learningrate = learningrate,
		momentum = momentum,
		learningrate_scale = learningrate_scale,
		hidden_dropout=hidden_dropout,visible_dropout=visible_dropout,
		output = output,
		W = W,
		vW = vW,
		B = B,
    vB = vB
	)
	
	m <- nrow(x);
	numbatches <- m / batchsize;
	s <- 0
	for(i in 1:numepochs){
		randperm <- sample(1:m,m)
		if(numbatches >= 1){
			for(l in 1 : numbatches){
			  s <- s + 1
				batch_x <- x[randperm[((l-1)*batchsize+1):(l*batchsize)], ] 
				if(is.vector(y)){
					batch_y <- y[randperm[((l-1)*batchsize+1):(l*batchsize)]] 
				}else if(is.matrix(y)){
					batch_y <- y[randperm[((l-1)*batchsize+1):(l*batchsize)], ] 
				}
				
				nn <- nn.ff(nn,batch_x,batch_y,s)
				nn <- nn.bp(nn)						
			}
		}
		#last fraction of sample
		if(numbatches > as.integer(numbatches)){
			batch_x <- x[randperm[(as.integer(numbatches)*batchsize):m], ]
			if(is.vector(y)){
				batch_y <- y[randperm[(as.integer(numbatches)*batchsize):m]]
			}else if(is.matrix(y)){
				batch_y <- y[randperm[(as.integer(numbatches)*batchsize):m], ]
			}
			s <- s + 1
			nn <- nn.ff(nn,batch_x,batch_y,s)
			nn <- nn.bp(nn)	
		}
			
		nn$learningrate <- nn$learningrate * nn$learningrate_scale;
	}
	
	nn
}

nn.ff <- function(nn,batch_x,batch_y,s){
  m <- nrow(batch_x)
  #do input dropout
  if(nn$visible_dropout > 0){
    nn$dropout_mask[[1]] <- dropout.mask(ncol(batch_x),nn$visible_dropout)
    batch_x <- t( t(batch_x) * nn$dropout_mask[[1]] )
  }
	nn$post[[1]] <- batch_x
	for(i in 2:(length(nn$size) - 1)){
	  nn$pre[[i]] <- t( nn$W[[i-1]] %*% t(nn$post[[(i-1)]])  + nn$B[[i-1]] )
		if(nn$activationfun == "sigm"){
			nn$post[[i]] <- sigm(nn$pre[[i]])
		}else if(nn$activationfun == "tanh"){
			nn$post[[i]] <- tanh(nn$pre[[i]])
		}else{
			stop("unsupport activation function!");
		}	
    if(nn$hidden_dropout > 0){
      nn$dropout_mask[[i]] <- dropout.mask(ncol(nn$post[[i]]),nn$hidden_dropout)
      nn$post[[i]] <- t( t(nn$post[[i]]) * nn$dropout_mask[[i]] )
    }
	}
	#output layer
	i <- length(nn$size)
  nn$pre[[i]] <- t( nn$W[[i-1]] %*% t(nn$post[[(i-1)]])  + nn$B[[i-1]] )
	if(nn$output == "sigm"){
		nn$post[[i]] <- sigm(nn$pre[[i]])
		nn$e <- batch_y - nn$post[[i]]
		nn$L[ s ] <- 0.5*sum(nn$e^2)/m
	}else if(nn$output == "linear"){
	  nn$post[[i]] <- nn$pre[[i]]
	  nn$e <- batch_y - nn$post[[i]]
	  nn$L[ s ] <- 0.5*sum(nn$e^2)/m
	}else if(nn$output == "softmax"){
	  nn$post[[i]] <- exp(nn$pre[[i]])
	  nn$post[[i]] <- nn$post[[i]] / rowSums(nn$post[[i]]) 
	  nn$e <- batch_y - nn$post[[i]]
	  nn$L[ s ] <- -sum(batch_y * log(nn$post[[i]]))/m
	}else{
	  stop("unsupport output function!");
	}	
  if(s %% 10000 == 0){
    message(sprintf("####loss on step %d is : %f",s,nn$L[ s ]))
  }
	
	nn
}


nn.bp <- function(nn){
	n <- length(nn$size)
  d <- list()
	if(nn$output == "sigm"){
		d[[n]] <- -nn$e * (nn$post[[n]] * (1 - nn$post[[n]]))
	}else if(nn$output == "linear" || nn$output == "softmax"){
		d[[n]] <- -nn$e
	}

	for( i in (n-1):2 ){
			if(nn$activationfun  == "sigm"){
				d_act <- nn$post[[i]] * (1-nn$post[[i]])
			}else if(nn$activationfun  == "tanh" ){
				d_act <- 1.7159 * 2/3 * (1 - 1/(1.7159)^2 * nn$post[[i]]^2)
			}
			d[[i]] <- (d[[i+1]] %*% nn$W[[i]]) * d_act
			if(nn$hidden_dropout > 0){
			  d[[i]] <- t( t(d[[i]]) * nn$dropout_mask[[i]] )
			}
	}
	
	for( i in 1:(n-1) ){
			dw <- t(d[[i+1]]) %*% nn$post[[i]] / nrow(d[[i+1]])
			dw <- dw * nn$learningrate
			if(nn$momentum > 0){
				nn$vW[[i]] <- nn$momentum * nn$vW[[i]] + dw
				dw <- nn$vW[[i]]
			}
			nn$W[[i]] <- nn$W[[i]] - dw
			
			db <- colMeans(d[[i+1]])
			db <- db * nn$learningrate
			if(nn$momentum > 0){
			  nn$vB[[i]] <- nn$momentum * nn$vB[[i]] + db
			  db <- nn$vB[[i]]
			}
			nn$B[[i]] <- nn$B[[i]] - db
	}
	nn
}

dropout.mask <- function(size,fraction){
  mask <- runif(size,min=0,max=1)
  mask[mask <= fraction] <- 0
  mask[mask > fraction] <- 1
  mask
}
