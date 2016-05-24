pcbic <-
function(eigenvals,n,pattern) {
     p <- length(eigenvals)
     l <- eigenvals
     k <- length(pattern)
     istart <- 1
     for(i in 1:k) {
          iend <- istart+pattern[i]
          l[istart:(iend-1)] = mean(l[istart:(iend-1)])
          istart <- iend
     }
     dev <- n*sum(log(l))
     dimen <- (p^2-sum(pattern^2))/2 + k
     bic <- dev + log(n)*dimen
     list(lambdaHat = l,Deviance = dev,Dimension = dimen,BIC = bic)
}
