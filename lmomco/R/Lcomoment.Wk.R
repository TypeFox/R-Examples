"Lcomoment.Wk" <-
function(k, r, n) {
   # Following notation of Serfling and Xiao (2006)
   #   compute the Wk weight factor for kth L-moment
   jn <- min(c(r-1,k-1))  # find the minimum for the loop end
   Wk <- sapply(0:jn, function(j) {
                          t2 <- lchoose(k-1,j); t3 <- lchoose(k-1+j,j)
                          t4 <- lchoose(n-1,j); t5 <- lchoose(r-1,j)
                          return( (-1)^(k-1-j) * (exp(t2 + t3 + t5 - t4)) )  })
   return(sum(Wk))
}
