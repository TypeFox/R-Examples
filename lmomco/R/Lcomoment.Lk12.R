"Lcomoment.Lk12" <-
function(X1, X2, k=1) {
  # Following notation of Serfling and Xiao (2006)
  #  compute the unbiased L-statistic estimator

  # Compute the concomitant of X2
  #    First sort X2 in ascending order, but need the indices
  #    Second rearrange X1 in the order of X2
  #I   <- sort(X2, decreasing=FALSE, index.return=TRUE)
  X12 <- X1[sort(X2, decreasing=FALSE, index.return=TRUE)$ix]
  n   <- length(X1)            # sample size
  SUM <- sum(sapply(1:n, function(r) { Lcomoment.Wk(k,r,n)*X12[r] } ))
  return(SUM/n)

  #sum <- 0                     # a summation
  #for(r in seq(1,n)) {         # for each value in the sample
  #  Wk  <- Lcomoment.Wk(k,r,n) # compute the weight factor
  #  sum <- sum + Wk*X12[r]     # sum them up
  #}                            # end of loop
  #Lk12 <- sum/n                # compute the expected value
  #return(Lk12)                 # return the L-comoment
}
