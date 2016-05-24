OptimClusts <-
function(P, Eps){
        ## P IS THE VECTOR CONTAINING THE AVERAGE SILHOUETTE WIDTHS
        ## Eps IS THE TUNING PARAMETER WHICH DEFINES MINIMUM DIFFERENCE BETWEEN MAXIMUM WIDTH AND THE NEXT LARGEST WIDTH. 
          N <- length(P);
          D <- seq(2,N+1, by = 1);
          MaxP <- max(P);
          S <- D[which(MaxP - P[-N] < Eps*MaxP)];
          Q = which.max(P)+1;
          S = S[S < Q];
          Z <- (MaxP - P[S-1])/(Q-S);
          if (length(S) > 0){
              return(S[which(Z == min(Z))]);
          }
          if (length(S) == 0){
              return(Q)
          }
}
