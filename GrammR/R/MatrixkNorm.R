MatrixkNorm <-
function(X,p){
        n = nrow(X);
        D <- matrix(0, ncol= n, nrow = n);
        for (i in 1:(n-1)){
                    for (j in (i+1):n){
                        if (is.finite(p)){
                                D[i,j] <- (sum(abs(X[i,] - X[j,])^p))^(1/p);
                        }
                        if (is.infinite(p)){
                                D[i,j] <- max(abs(X[i,] - X[j,]));
                        }
                    D[j,i] <- D[i,j];
                  }
        }
        return(D);
}
