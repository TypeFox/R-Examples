    #second component
eblup.mse.f.c2 <-
    function(gamma.i
             , X.i#design matrix for sampled elements in domain
             , X.bar.i #mean of pop. design matrix
             , sum.A.i #sum of all m A.i matrices
             ,...){
            #7.2.12
                #mean of the design matrix of the sampled elemets
        x.bar.i <- as.matrix(apply(X.i, 2, mean))
        x.tmp <- (X.bar.i - gamma.i * x.bar.i) #%*%
        res <- t(x.tmp) %*% solve(sum.A.i) %*% x.tmp
        return(res)
    }

