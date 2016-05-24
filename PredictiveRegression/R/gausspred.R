# The Gauss predictor.
# Description: paper by Vovk et al. in the Annals of Statistics, 2009
# (full version: arXiv ST/0511522), referred to as "paper".
# This is the standard textbook predictor.

gausspred <- function(train,test,epsilons=c(0.05,0.01)) {

# train is the matrix of training observations
#   each row is the vector of explanatory variables followed by the response
# test is the matrix of explanatory variables for test observations
#   each row is the vector of explanatory vatiables
# epsilons is the list of significance levels to consider

  N <- dim(train)[1];   # number of training observations
  K <- dim(train)[2]-1; # number of explanatory variables in the training set
  N2 <- dim(test)[1];   # number of test observations
  K2 <- dim(test)[2];   # number of explanatory variables in the test set
  flag <- 0;            # normal termination of the program

  up <- array(Inf,c(N2,length(epsilons)));    # initializing upper limits
  low <- array(-Inf,c(N2,length(epsilons)));  # initializing lower limits

  if (K2 != K) {  # illegal parameters
    flag <- 1;    # code 1 (illegal parameters)
    return(list(low,up));
  }

  if (N < K+3) {    # in this case computation is impossible (Table 1 in the paper)
    flag <- 2;      # code 2 (too few observations)
    return(list(low,up));
  }

  # training set:
  ZZ <- train[,1:K];   # matrix of explanatory variables
  dim(ZZ) <- c(N,K);
  ZZ <- cbind(1,ZZ);   # extended matrix of explanatory variables
  dim(ZZ) <- c(N,K+1);
  yy <- train[,K+1];   # vector of responses
  dim(yy) <- c(N,1);

  # test set:
  ZZ2 <- test[,1:K];
  dim(ZZ2) <- c(N2,K);
  ZZ2 <- cbind(1,ZZ2);  # extended matrix of explanatory variables
  dim(ZZ2) <- c(N2,K+1);

  # computing the prediction intervals: start
  toinvert <- t(ZZ)%*%ZZ;                   # the matrix to invert (not explicitly)
  gamma_hat <- solve(toinvert,t(ZZ)%*%yy);  # the maximum likelihood weights
  center <- ZZ2%*%gamma_hat;                # the point least-squares predictions

  sigma_hat_sq <- sum((yy-ZZ%*%gamma_hat)^2) / (N-K-1); # estimated variance
  VV <- array(0,dim=c(N2,length(epsilons)));  # the prediction interval semi-widths initialized

  for (n in 1:N2)     # loop over the test set: start
  {
    zz <- ZZ2[n,];    # the explanatory variables of the new observation
    dim(zz) <- c(K+1,1);
    # the semi-width before taking epsilon into account:
    VV_before <- sqrt((1+t(zz)%*%solve(toinvert,zz)) * sigma_hat_sq);
    for (eps_index in 1:length(epsilons))  # the epsilons loop: start
    {
      epsilon <- epsilons[eps_index];
      VV <- qt(1-epsilon/2,N-K-1) * VV_before;  # prediction interval semi-width
      low[n,eps_index] <- center[n]-VV;
      up[n,eps_index] <- center[n]+VV;
    }                                      # the epsilons loop: end
  }                   # loop over the test set: end
  list(low,up,flag);
}
