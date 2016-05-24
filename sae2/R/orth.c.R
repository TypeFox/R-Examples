# function to compute an orthogonal complement of X, used in computing
# the constrained likelihood

orth.c <- function(X,tol=.1e-7) {
     nr <- nrow(X)
     nc <- ncol(X)
     XE <- cbind(X,diag(nr))
     use.c <- rep(TRUE,times=nr+nc)
     for (i in 1:(nr+nc)) {
       w<- sum(XE[,i]^2)
       if (w < tol) use.c[i] <- FALSE else
         {XE[,i] <- XE[,i]/sqrt(w)
          if (i < nr+nc) {
            for (j in (i+1):(nr+nc)) {
                XE[,j] <- XE[,j] - sum(XE[,i]*XE[,j]) * 
                          XE[,i]
            }
          }   
         }
       }
      cols <- (nc+1):(nc+nr)
      cols <- cols[use.c[cols]==TRUE]
      XE <- XE[,cols] 
      return(XE)
      }
