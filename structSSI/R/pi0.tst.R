  # This function computes the estimate of pi0 for
  # each group using the TWo-Stage method proposed
  # by Benjamini, Krieger, and Yekutieli (2006) and
  # explained by Xu, Zhou, and Zhao (2010). The idea
  # is that it estimates the number of alternatives
  # in a group of p-values by the numberof rejections
  # made by the Classical Benjamini-Hochberg Procedure.
  #
  # Input: 1) p.val - [numeric vector] - A vector of the unadjusted p-values
  # associated with a particular group in the Group Benjamini-Hochberg
  # procedure.
  # 2) alpha - [numeric] - The level of the overall multiple testing
  # procedure. In our case, the size of the GBH procedur
  # that we are working with. This is just a number.
  #
  # Output: 1) pi0 - A numeric estimate of the proportion
  # of null hypotheses within the group of p-values
  # provided. This is just some fraction.
  # Note: This function uses the Benjamini-Hochberg procedure,
  # so it requires multtest.

  pi0.tst <- function(p.val, alpha = 0.05){
    alpha.prime <- alpha/(1 + alpha)
    n_g <- length(p.val)
    adjustment <- mt.rawp2adjp(p.val, proc = "BH")
    rejected <- mt.reject(adjustment$adjp, alpha.prime)
    n.rejected <- rejected$r[,2]
    estimate <- (n_g - n.rejected)/n_g
    return(estimate)
  }
