mvlognormal <-
function(n = 1,Mu,Sigma, R){
      ## Generate Lognormal random variables with zeros. 
      ## Mu - Mean of the actual lognormal variables
      ## Sigma - Diagonal of the covariance matrix of the actual variables.
      ## R - Correlation matrix for the log-transformed normal variables. This is to ensure the normal distribution corresponding to the lognormals exists
      
      if (dim(R)[1] != dim(R)[2]){
	  	stop("Correlation matrix is not square.");
      }
      p <- length(Mu);
      if (length(Sigma) != p){
	  	stop("Mean and covariance matrix are not of same dimension")
      }
      Alpha <- matrix(0,p,1);
      Beta <- matrix(0,p,p);
      
      ## Alpha and Beta are the converted mean and covariance-diagonal matrix of the log-transformed normal variable.
      for (i in 1:p){
		  if (abs(Mu[i]) >= .Machine$double.eps){
			    Alpha[i] = log(Mu[i]) - (1/2)*log(1 + Sigma[i]/Mu[i]^2);
			    Beta[i,i] = log(1 + Sigma[i]/(Mu[i]^2));
		  }
		  if (abs(Mu[i]) < .Machine$double.eps){
			    Alpha[i] = -Inf;
			    Beta[i,i] = 0;
		  }
      }
      Delta = sqrt(Beta)%*%R%*%sqrt(Beta);
      Delta.Eigen = eigen(Delta);
      Delta.EigenVal <- Delta.Eigen$values*(abs(Delta.Eigen$values) >= 1e-12);
      RootDelta = (Delta.Eigen$vectors)%*%diag(sqrt(Delta.EigenVal))%*%t(Delta.Eigen$vectors);
      X <- matrix(rnorm(p*n,0,1), nrow = n, ncol = p);
      for (i in 1:n){
	  	X[i,] <- exp(Alpha + RootDelta%*%X[i,]);
      }      
      return(X);
}
