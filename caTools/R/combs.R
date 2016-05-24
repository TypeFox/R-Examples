#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

combs = function(v,k) {
# combs(V,K) - finds all unordered combinations of K elements from vector V 
#  V is a vector of length N
#  K is a integer 
# combs(V,K) creates a matrix with N!/((N-K)! K!) rows
# and K columns containing all possible combinations of N elements taken K at a time.
# example: combs(1:3,2) returns matrix with following rows (1 2), (1 3), (2 3)
  n = length(v)
  if      (n==k  ) P = matrix(v,1,n)
  else if (k==1  ) P = matrix(v,n,1)
  else if (k==n-1) P = matrix( rep(v, each=n-1), n, n-1)
  else if (k< n) {
    P = matrix(0,0,k)
    if (k < n & k > 1) {
      for (i in 1:(n-k+1)) {
        Q = combs(v[(i+1):n],k-1)
        j = nrow(Q)
        P = rbind(P, cbind(rep(v[i],j), Q))
      }
    }
  } else 
    stop("combs: number m has to be smaller or equal to length of vector v")
  return(P)
}

