
  # This method implements the Least Slope Estimator proposed by
  # Benjamini-Hochberg (2000) and described by Xu, Zhou, and Zhao
  # in their Group Benjamini-Hochberg paper.
  #
  # 1) For p-values ingroup g, starting from i = 1, compute
  # \code{l\_{g},i = (n\_{g} + 1 - i)/(1 - P\_{g,i})} where n\_g is the number
  # of p-values in group g and P_g,i is the ith ordered p-value
  # in group g.
  # 2) As i increases, stop when \code{l\_{g},j > l\_{g,(j - 1)}} for the first
  # time.
  # 3) For each group, compute the LSL estimator of pi_{g,0} as
  # \code{gamma\_{g, lsl} = min((floor(l\_{g},j) + 1)/n\_{g, 1}).
  #
  # Here are details of this function.
  # Input: 1) p.val [numeric vector] - A vector of p-values.
  # Output: 1) pi0 - An estimate for the proportion of null
  # hypotheses within the set of p-values provided.

pi0.lsl <- function(p.val){
    p.val <- sort(p.val)
    n_g <- length(p.val)

    i <- 1
    while(TRUE){
        if(i >= 2){
            l_g.i.prev <- l_g.i
        } else {
            l_g.i.prev <- 10000   # We don't want to stop on the first iteration, so
                                        # we set the values used to estimate pi0 very high
        }
        if(p.val[i] < 1){
            l_g.i <- (n_g + 1 - i)/(1 - p.val[i])
        } else {
            return(l_g.i.prev) # If you're dividing by zero, then you guarantee an increase
        }
        if(l_g.i > l_g.i.prev || i == length(p.val)){
            pi0 <- (floor(l_g.i) + 1)/n_g
            pi0 <- min(pi0, 1)
            return(pi0)
        }
        i <- i + 1
    }
}
