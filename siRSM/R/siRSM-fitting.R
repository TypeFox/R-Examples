#################################################################################
# single.run: where alternating optimization is implemented 
#################################################################################
single.run <- function (data, w0, int.only=FALSE, eps=0.005, precision=1e-6, max.iter=10000, debug=FALSE, trace=TRUE)
{
# Written by Huan Cheng, 2014
# Modified by Mu Zhu, 2014
# data: must be structured so that 
#   1st column = y; 
#   next K columns are candidates for factor one;
#   next K columns are candidates for factor two;
# eps: stepsize;
# precision: stopping criterion;
# max.iter: max number of iterations;

  #### extract information from the dataset ####
  
  # N: the sample size
  N <- nrow(data)  
  # K: dimensionality of w
  K <- (ncol(data)-1)/2 
  name <- names(data)
  # X = candidates for factor one
  # Z = candidates for factor two
  # y = response
  X <- as.matrix(data[,2:(K+1)]) 
  Z <- as.matrix(data[,(K+2):(2*K+1)])
  y <- as.vector(data[,1])
  
  if(missing(w0)) 
  {
    w0 <- runif(K,-1,1) 
  }
  # standardize
  w_new <- w0/sqrt(sum(w0^2)) 
  
  #### inialize the paramters ####
  if (int.only){
    b_new <- b <- runif(4)
  }
  else{
    b_new <- b <- runif(6)
  }
  # the difference of values between two steps
  diff_b <- diff <- 1
  # number of iteration
  count <- 0 
  loss <- c()
  dx <- matrix(0,ncol=N,nrow=K)
  dmz=matrix(0,N,K)
  
  #### gradient descent iterations ####
  
  while((count< max.iter) && (diff > precision) && (diff_b > precision))      
  {
    w <- w_new
    u <- X%*%w # N by 1
    v <- Z%*%w # N by 1 
  
    ### estimate the response surface ###
    
	if (int.only){
	  b_new <- coef(lm(y~u+v+I(u*v)))
	}
	else{
      b_new <- coef(lm(y~u+v+I(u^2)+I(u*v)+I(v^2)))
	}
    diff_b <- sum(abs(b_new-b))
    b <- b_new 
    
    ### update w ###
    
 	if (int.only){
      x.design <- cbind(rep(1,N),u,v,u*v)         # N by 4
	}
 	else{
      x.design <- cbind(rep(1,N),u,v,u^2,u*v,v^2) # N by 6 
	}
    res <- y - x.design%*%b 
    

    for (i in 1:K)
    {
	  if (int.only) {
      dx[i,] <- t(-b)%*%t(cbind(rep(0,N),
                                X[,i],
                                Z[,i],
                                (X%*%w*Z[,i]+Z%*%w*X[,i])))      
	  }
	  else{
      dx[i,] <- t(-b)%*%t(cbind(rep(0,N),
                                X[,i],
                                Z[,i],
                                (2*X%*%w*X[,i]),
                                (X%*%w*Z[,i]+Z%*%w*X[,i]),
                                (2*Z%*%w*Z[,i])))      
	  }
    }
    # take derivative of MSE with respect to w 
    dL<- 2* dx %*% res/N 

	if (debug) { 	# MZ verification code ... used while debugging
	                # does not apply to case of int.only=TRUE
    cat('Huan gradient\t', as.vector(dL), '\n')	
	for (i in 1:N)
	{
 	  dmz[i,]=t(b)%*% rbind(
	      rep(0,K),
	      X[i,],
		  Z[i,],
		  2*t(w)%*%X[i,]%*%t(X[i,]),
		  t(w)%*%X[i,]%*%t(Z[i,]) + t(w)%*%Z[i,]%*%t(X[i,]),
		  2*t(w)%*%Z[i,]%*%t(Z[i,]))
	}	  
	dL.check=(-2)*t(dmz) %*% res/N
    cat('Mu gradient\t', as.vector(dL.check), '\n')
	
	J=dmz
	dGN=(-1)*solve(t(J)%*%J)%*%t(dmz) %*% res/N
    cat('Gauss-Newton direction\t', as.vector(dGN), '\n')
    }
	
    # update w                       
    w_new <- w-(eps*dL) 
    w_new <- w_new/((sum(w_new^2))^0.5) 
    if(w_new[1]<0) w_new <- -w_new 
	if (trace) {
		cat(count, ': w\t', as.vector(w_new), '\n')
	}
    diff<- sum(abs(w-w_new))
    loss[count+1]<- sum(res^2)/N
    count<- count+1 
  }
  
  if (int.only){
    reg <- lm(y~u+v+I(u*v))
  }
  else{
    reg <- lm(y~u+v+I(u^2)+I(u*v)+I(v^2))
  }
  reg$call <- summary(reg)
  
  return<- list(w=w, 
                coefficients=b,
				u=as.vector(u), v=as.vector(v), y=as.vector(y),
				type=ifelse(int.only,"interaction.only","full.quadratic"),
                iter=count,
                residuals=res,
                mse.path=loss,
                mse=loss[count-1],
                model=reg$call)  
}

#################################################################################
# multi.run: calls single.run with different initial values 
#################################################################################
multi.run <- function(y, X, Z, rep, interaction.only=FALSE, use.parallel=TRUE)
{ 
  if (ncol(X) != ncol(Z)) {
    cat('Covariate dimensions do not match!\n')
  }
  else{
    data=cbind(y, X, Z)
    K <- ncol(X)
    # rep is the number of intial values tried, default=2*dimension set by MZ
    if (missing(rep)) rep=2*K

   if (use.parallel) {   # new parallel code using foreach, etc
    fit <- as.list(numeric(rep))
    loss <- numeric(rep)
	MAX=detectCores()-1
    cl <- makeCluster(ifelse(rep>MAX,MAX,rep),type='SOCK') 
	cat('Using parallel processing ...\n')
	print(cl)
	cat('\n')
    registerDoSNOW(cl=cl)
	clusterExport(cl, 'single.run')
     rslt <- foreach(i=1:rep, .verbose=FALSE) %dopar% {
      w0 <- runif(K,-1,1);
      f <- single.run(data, w0, int.only=interaction.only, trace=FALSE);
	  return(list(fit=f, loss=f$mse));
    }
	stopCluster(cl)
	for (i in 1:rep){
	  fit[[i]]=rslt[[i]]$fit
	  loss[[i]]=rslt[[i]]$loss
	}
   }
   else{     # old serial code
    w0 <- matrix(0,nrow=rep,ncol=K)
    fit <- as.list(numeric(rep))
    loss <- numeric(rep)
    for(i in 1:rep)
    {
      cat("trying different initial values, trial", i, "... ")
      # uniformly generate intial values for 2 
      w0[i,] <- runif(K,-1,1)
      # call function single.run() 
      f <- single.run(data, w0[i,],trace=FALSE)
      fit[[i]] <- f
      loss[i] <- f$mse   
      cat("MSE = ",round(loss[i],4),"\n")
    }
   }

    # choose the one with minimum mse 
    po <- which.min(loss)
    call <- fit[[po]]
    return(call)
  }
}
siRSM.default <- function(y, U, V, trial, interaction.only=FALSE, use.parallel=TRUE)
{
  # original code uses y as response, X, Z as two sets of covariates
  # changed (X,Z) to (U,V) in the main call so that potential users in psychology
  # won't get confused, but didn't change in main call to reduce
  # mistakes
  tmp <- multi.run(y=y, X=U, Z=V, rep=trial, interaction.only=interaction.only, use.parallel=use.parallel)  
  class(tmp) <- "siRSM"
  tmp
}
siRSM <- function(y, U, V, trial, interaction.only=FALSE, use.parallel=TRUE) UseMethod("siRSM")

print.siRSM.default <- function(x, ...)
{
  cat('Single Index Response Surface Model\n\n')
  cat('***Single Index Estimate***\n\n', round(as.vector(x$w),4), '\n\n')
  cat('***Corresponding Quadratic Surface***\n')
  print(x$model)  
}
print.siRSM <- function(x, ...) UseMethod("print.siRSM")
