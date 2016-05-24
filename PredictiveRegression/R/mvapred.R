# The MVA predictor.
# Description: paper by Vovk et al. in the Annals of Statistics, 2009
# (full version: arXiv ST/0511522), referred to as "paper".

mvapred <- function(train,test,epsilons=c(0.05,0.01),ridge=0) {

# train is the matrix of training observations
#   each row is the vector of explanatory vatiables followed by the response
# test is the matrix of explanatory variables for test observations
#   each row is the vector of explanatory vatiables
# ridge is the ridge coefficient (usually a small positive number, such as 0.01)
# epsilons is the list of significance levels to consider

  N <- dim(train)[1];   # number of training observations
  K <- dim(train)[2]-1; # number of explanatory variables in the training set
  N2 <- dim(test)[1];   # number of test observations
  K2 <- dim(test)[2];   # number of explanatory variables in the test set
  flag <- 0;            # normal termination of the program

  up <- array(Inf,c(N2,length(epsilons)));   # initializing upper limits
  low <- array(-Inf,c(N2,length(epsilons))); # initializing lower limits

  if (K2 != K) {  # illegal parameters
    flag <- 1;    # code 1 (illegal parameters)
    return(list(low,up,flag));
  }

  if (N < 3) {  # in this case computation is impossible (Table 1 in the paper)
    flag <- 2   # code 2 (too few observations)
    return(list(low,up));
  }

  # training set:
  ZZ <- train[,1:K];    # matrix of explanatory variables
  dim(ZZ) <- c(N,K);
  ZZ <- cbind(1,ZZ);    # extended matrix of explanatory variables
  dim(ZZ) <- c(N,K+1);
  yy <- train[,K+1];    # vector of responses
  dim(yy) <- c(N,1);
  # test set:
  ZZ2 <- test[,1:K];
  dim(ZZ2) <- c(N2,K);
  ZZ2 <- cbind(1,ZZ2);  # extended matrix of explanatory variables
  dim(ZZ2) <- c(N2,K+1);

  for (n in 1:N2)     # The loop over the test set: start
  {
    ZZ_ext <- rbind(ZZ,ZZ2[n,]);
    dim(ZZ_ext) <- c(N+1,K+1);
    yy_ext <- rbind(yy,0);
    dim(yy_ext) <- c(N+1,1);
    OO1 <- rbind(array(0,c(N,1)),1);
    dim(OO1) <- c(N+1,1);
    toinvert <- t(ZZ_ext) %*% ZZ_ext + ridge*diag(K+1);
    A <- yy_ext - ZZ_ext %*% solve(toinvert,t(ZZ_ext)) %*% yy_ext;
    B <- OO1 - ZZ_ext %*% solve(toinvert,ZZ_ext[N+1,]);
    A <- A - mean(A[1:N]);
    B <- B - mean(B[1:N]);
    AA <- sum(A[1:N]^2);
    BB <- sum(B[1:N]^2);
    AB <- sum(A[1:N]*B[1:N]);

    # Solving the quadratic inequality ay^2 + 2by + c < 0;
    # solutions: -b/a \pm \sqrt(D)/a, where D = b^2 - ac.

    for (eps_index in 1:length(epsilons))   # the epsilons loop: start
    {
      epsilon <- epsilons[eps_index];
      t_sq <- qt(1-epsilon/2,N-1)^2;   # quantile squared
      # coefficients of the quadratic equation
      a <- N * (N-1) * B[N+1]^2 - t_sq * (N+1) * BB;
      b <- N * (N-1) * A[N+1]*B[N+1] - t_sq * (N+1) * AB;
      c <- N * (N-1) * A[N+1]^2 - t_sq * (N+1) * AA;

      lower <- -Inf;
      upper <- Inf;

      if (a > 0){
        center <- -b/a;    # center of the prediction interval
        D <- b^2 - a*c;    # discriminant
        if (D >= 0) {
          D_sqrt <- sqrt(D);
          upper <- center + D_sqrt/a;
          lower <- center - D_sqrt/a;
        } else {   # empty prediction
          upper <- -Inf;
          lower <- Inf;
        }
      } else if (a < 0) {  # prediction = union of two rays
        upper <- Inf;
        lower <- -Inf;
      } else {             # a=0, and so ay^2 + 2by + c is a straight line
        boundary <- -c/(2*b); # the boundary of the prediction interval (which is a ray if b!=0)
        if (b > 0) {
          upper <- boundary;
          lower <- -Inf;
        } else if (b < 0) {
          upper <- Inf;
          lower <- boundary;
        } else {    # a=0 and b=0
          if (c<0) { # prediction interval = real line
            upper <- Inf;
            lower <- -Inf;
          } else {  # empty prediction
            upper <- -Inf;
            lower <- Inf;
          }
        }
      }

      up[n,eps_index] <- upper;
      low[n,eps_index] <- lower;
    }  # the epsilons loop:end
  }  # the loop over the test: end
  list(low,up,flag);
}
