"binomial.expand" <-
function(Y,k,w){
      
    if (is.matrix(Y)){
            r <- rep(Y[,1],k)
            n <- rep(Y[,2]+Y[,1],k)
    } else { N <- NROW(Y)
             r <- rep(as.numeric(Y*w),k)
             n <- rep(w,k)
             #w<-rep(1,N)
    }
    Y <- cbind(r,"n-r"=(n-r))
    PY<-Y[,1]/(Y[,1]+Y[,2])     
    list(Y,PY,r,n)
}

