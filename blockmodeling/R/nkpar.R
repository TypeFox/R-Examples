"nkpar" <-
function(n, k) {
# Author: Chris Andrews <candrews@buffalo.edu>
     sum((-1)^seq(0,k-1) * choose(k, seq(0,k-1)) * (k-seq(0,k-1))^n) / factorial(k)
}

