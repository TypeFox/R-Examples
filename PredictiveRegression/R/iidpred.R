# The IID predictor
# (= Ridge Regression Confidence Machine
# simplified to compute prediction intervals only).
# Description: paper by Vovk et al. in the Annals of Statistics, 2009
# (full version: arXiv ST/0511522), referred to as "paper".
# The algorithm implemented by the IID predictor
# is also described in: Vovk, Gammerman, Shafer
# "Algorithmic Learning in a Random World", Springer, New York, 2005,
# pages 30 - 34.

iidpred <- function(train,test,epsilons=c(0.05,0.01),ridge=0) {

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

  if (K2 != K) {   # illegal parameters
    flag <- 1;     # code 1 (illegal parameters)
    return(list(low,up,flag));
  }

  # training set:
  ZZ <- train[,1:K];    # matrix of explanatory variables
  dim(ZZ) <- c(N,K);
  ZZ <- cbind(1,ZZ);    # extended matrix of explanatory variables (paper, (4))
  dim(ZZ) <- c(N,K+1);
  yy <- train[,K+1];    # vector of responses
  dim(yy) <- c(N,1);
  # test set:
  ZZ2 <- cbind(1,test); # extended matrix of explanatory variables (paper, (4))
  dim(ZZ2) <- c(N2,K+1);

  A <- array(0,c(N+1,1));           # vector a from the IID predictor (the paper, Appendix B)
  B <- array(0,c(N+1,1));           # vector b from the IID predictor
  toinvert <- array(0,c(K+1,K+1));  # auxiliary matrix
  inverseZZ <- array(0,c(K+1,N+1)); # auxiliary matrix

  P <- array(0,c(2*N+2,1));   # array P of critical points
  NM <- array(0,c(2*N+2,1));  # combination of the arrays N and M from RRCM
                              # its size is sizeP

  for (n2 in 1:N2) {    # the loop over the test set: start
    P[1] <- -Inf;   # the first critical point is -infty
    sizeP <- 1;     # size of the array P: initialized

    ZZ_ext <- rbind(ZZ,ZZ2[n2,]);  # training set extended by adding the test observation
    dim(ZZ_ext) <- c(N+1,K+1);
    yy0 <- rbind(yy,0);            # response vector extended by adding 0
    dim(yy0) <- c(N+1,1);
    OO1 <- rbind(array(0,c(N,1)),1); # vector 0...01
    dim(OO1) <- c(N+1,1);

    toinvert <- t(ZZ_ext)%*%ZZ_ext + ridge*diag(K+1);
    inverseZZ <- solve(toinvert,t(ZZ_ext));   # <- inverse %*% t(ZZ_ext)
    A <- yy0 - ZZ_ext %*% (inverseZZ %*% yy0);
    B <- OO1 - ZZ_ext %*% inverseZZ[,N+1];

    # make sure B >= 0:
    A[B<0] <- -A[B<0];
    B[B<0] <- -B[B<0];
    # initialize the counts of the number of S_n, n=1,...,N, covering -infty and infty:
    L <- 0;
    R <- 0;

    # finding the critical points and NM: start
    for (n in 1:N) {
      if (B[n] != B[N+1]) {   # 2 critical points in this case (which may coincide)
        point1 <- (A[n]-A[N+1]) / (B[N+1]-B[n]);
        point2 <- -(A[n]+A[N+1]) / (B[N+1]+B[n]);
        P[sizeP+1] <- min(point1,point2);
        P[sizeP+2] <- max(point1,point2);
        if (B[n] < B[N+1]) {  # the case of interval
          NM[sizeP+1] <- 1;
          NM[sizeP+2] <- -1;
        } else {     # B[n] > B[N+1]: the case of two rays
        # cat("Two rays\n");  # testing (this point is sometimes reached for the data set used in the paper)
          NM[sizeP+1] <- -1;
          NM[sizeP+2] <- 1;
          L <- L+1;
          R <- R+1;
          if (point1 == point2) { # the rays are not disjoint
            sizeP <- sizeP-2;     # reverse adding the critical points
          }
        }
        sizeP <- sizeP+2;
      } else {     # B[n]==B[N+1]
        # cat("Equal slopes\n");  # testing (this point is never reached for the data set used in the paper)
        if (A[n] == A[N+1]) {     # the trivial case
          L <- L+1;
          R <- R+1;
        } else {   # B[n]==B[N+1] and A[n]!=A[N+1]
          if (B[N+1] != 0) {  # the case of ray; 1 critical point in this case
            point1 <- -(A[n]+A[N+1]) / (B[N+1]+B[n]);
            if (A[n] > A[N+1]) {  # the ray is going to the right
              NM[sizeP+1] <- 1;
              R <- R+1;
            } else {              # the ray is going to the left
              NM[sizeP+1] <- -1;
              L <- L+1;
            }
            sizeP <- sizeP+1;
          } else {                # the lines a+by are horizontal
            if (abs(A[n]) >= abs(A[N+1])) {
              L <- L+1;
              R <- R+1;
            }
          }
        }
      }
    }
    P[sizeP+1] <- Inf;    # the last critical point
    sizeP <- sizeP+1;
    NM[1] <- L+1;
    NM[sizeP] <- -R-1;
    # finding the critical points and NM: end
    # cat("Should be zero: ",sum(NM),"\n");   # testing (sum(NM) is always zero)

    # ordering P and NM:
    P_order <- order(P,-NM);
    P <- P[P_order];
    NM <- NM[P_order];

    eps_order <- order(epsilons);

    p <- 0;   # current p-value
    P_reached <- 1;
    for (eps_index in 1:length(epsilons)) {  # the epsilons loop for the lower bounds: start
      eps_real_index <- eps_order[eps_index];
      epsilon <- epsilons[eps_real_index];
      found <- FALSE;
      for (P_index in P_reached:sizeP) {
        p <- p + NM[P_index]/(N+1);
        if (p > epsilon) {
          found <- TRUE;
          low[n2,eps_real_index] <- P[P_index];
          P_reached <- P_index;
          p <- p - NM[P_index]/(N+1);   # preventing double counting
          break;
        }
      }
      if (!found) {
        for (eps_ind in eps_index:length(epsilons))	{
          eps_real_ind <- eps_order[eps_ind];
          low[n2,eps_real_ind] <- Inf;
          up[n2,eps_real_ind] <- -Inf;
        }
        break;
      }
    }                                       # the epsilons loop for the lower bounds: end

    p <- 0;     # current p-value
    P_reached <- sizeP;
    for (eps_index in 1:length(epsilons)) { # the epsilons loop for the upper bounds: start
      eps_real_index <- eps_order[eps_index];
      epsilon <- epsilons[eps_real_index];
      found <- FALSE;
      for (P_index in P_reached:1) {
        p <- p - NM[P_index]/(N+1);
        if (p > epsilon) {
          found <- TRUE;
          up[n2,eps_real_index] <- P[P_index];
          P_reached <- P_index;
          p <- p + NM[P_index]/(N+1);  # preventing double counting
          break;
        }
      }
      if (!found) {
        cat("You have found a bug in the program.\n");
        cat("Please contact the package's maintainer.\n");
      }
    }                                    # the epsilons loop for the upper bounds: end
  }    # the loop over the test set: end

  # setting the flag:
  max_eps <- max(epsilons);		# the largest significance level
  if (N+1 < 1/max_eps) flag <- 2;	# too few observations for all significance levels

  list(low,up,flag);  # output of the function
}
