covariance <-
function(x,stat,method,...){
p <- ncol(x) # quality characteristics
m <- nrow(x) # sample
if (class(x) == "matrix" || class(x) == "data.frame") (x <- array(data.matrix(x),c(m,p,1)))
n <- dim(x)[3] # observations or sample size

s.jk <- matrix(0,m,p ^ 2) # matrix of the covariances of observations
SS <- matrix(0,m,1) # matrix of /S/ statistic 

if(n > 1){
     arrays <- expand.grid(1:p,1:p)
 
     for (i in 1 : m){
         for(j in 1 : p ^ 2){
     s.jk[i,j] <- cov(x[i,arrays[j,1],],x[i,arrays[j,2],])
        }
    } 

     S <- matrix(colMeans(s.jk),p,p)

     for (ii in 1 : m){
         SS[ii] <- det(matrix(s.jk[ii,],p,p))
        }

     if(missing(stat)) (return(S))
     else (return(SS))

    }    

if(n == 1){
 if(missing(method))(method="sw")
 
 if(method == "sw"){
     B <- matrix(0,p,p)
     w <- sweep(x,2,(apply(x,2,mean))) #compute de value minus the mean
 for(i in 1:m){
     B <- B + w[i,,] %*% t(w[i,,])
    }
 S <- s1 <- B/(m - 1)
 }
  
 if(method == "hm"){
     V <- matrix(0,m-1,p)
     for(i in 1:m-1){
     V[i,] <- x[i+1,,] - x[i,,]
    }
 S <- s2 <- .5 * t(V) %*% V / (m - 1)
    }
 

 return(S)
}


}
