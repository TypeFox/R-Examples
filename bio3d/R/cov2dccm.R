cov2dccm <- function(vcov, method = c("pearson", "lmi"), ncore = NULL) {
   method = match.arg(method)

   ncore = setup.ncore(ncore)
   if(ncore == 1) mclapply = lapply
   
   x <- vcov
   ccmat = switch(method, 
      pearson = {
          n <- nrow(x)
          np <- pairwise(n/3)
       
          d <- sqrt(colSums(matrix(diag(x), nrow=3)))
       
          ltv <- mclapply(1:nrow(np), function(i) {
             i1 <- (np[i, 2] - 1) * 3 + 1
             i2 <- (np[i, 1] - 1) * 3 + 1
             sum(diag(x[i1:(i1+2), i2:(i2+2)]))/ # sum of diagnol of submatrix
                (d[np[i, 2]] * d[np[i, 1]])   # divided by product of standard deviations
          } )
       
          ccmat <- matrix(0, n/3, n/3)
          ccmat[lower.tri(ccmat)] <- unlist(ltv)
       
          # make full matrix
          ccmat <- ccmat + t(ccmat)
          diag(ccmat) <- 1
          ccmat
      },
      lmi = {
      # rm:r-value matrix
          cm <- x
          l <- dim(cm)[1]/3
          rm <- matrix(nrow=l, ncol=l)
          d <- 3
          ij <- pairwise(l)
      
      # list1: marginal-covariance 
          list1 <- mclapply(1:l, function(i) det(cm[(3*i-2):(3*i), (3*i-2):(3*i)]) )
          dm <- unlist(list1)
      
      # list2: pair-covariance
          list2 <- mclapply(1:nrow(ij), function(i) {
              x <- det(cm[c((3*ij[i,1]-2):(3*ij[i,1]),(3*ij[i,2]-2):(3*ij[i,2])), c((3*ij[i,1]-2):(3*ij[i,1]),(3*ij[i,2]-2):(3*ij[i,2]))])
              y <- 1/2 * (log(dm[ij[i,1]]) + log(dm[ij[i,2]]) - log(x))
              (1 - exp(-2 * y / d))^(1/2)
              }
          )
          list2 <- unlist(list2)
      
          for (k in 1:nrow(ij)) {
              rm[ij[k, 1], ij[k, 2]] <- list2[k]
          }
          rm[lower.tri(rm)] = t(rm)[lower.tri(rm)]
          diag(rm) <- 1
          rm
      } )
   class(ccmat)=c("dccm","matrix")
   return(ccmat)
}
