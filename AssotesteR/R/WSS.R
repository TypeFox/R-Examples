#'WSS: Weighted Sum Statistic
#'
#'The WSS method has been proposed by Madsen and Browning (2009) as a pooling
#'approach. In WSS, rare-variant counts within the same gene for each
#'individual are accumulated rather than collapsing on them. Second, it
#'introduces a weighting term to emphasize alleles with a low frequency in
#'controls. Finally, the scores for all samples are ordered, and the WSS is
#'computed as the sum of ranks for cases. The significance is determined by a
#'permutation procedure.
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem wss.stat wss statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{CMC}}
#'@references Madsen BE, Browning SR (2009) A Groupwise Association Test for
#'Rare Mutations Using a Weighted Sum Statistic. \emph{PLoS Genetics},
#'\bold{5(2)}: e1000384
#'@examples
#'
#'  \dontrun{
#'  
#'  # number of cases
#'  cases = 500
#'
#'  # number of controls
#'  controls = 500
#'
#'  # total (cases + controls)
#'  total = cases + controls
#'
#'  # phenotype vector
#'  phenotype = c(rep(1, cases), rep(0, controls))
#'
#'  # genotype matrix with 10 variants (random data)
#'  set.seed(123)
#'  genotype = matrix(rbinom(total*10, 2, 0.05), nrow=total, ncol=10)
#'
#'  # apply WSS with 500 permutations
#'  mywss = WSS(phenotype, genotype, perm=500)
#'  mywss
#'  }
#'
WSS <-
function(y, X, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  
  ## running wss method
  wss.stat = my_wss_method(y, X)  
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_wss_method(y[perm.sample], X) 
    }
    # p-value 
    perm.pval = sum(x.perm >= wss.stat) / perm
  }
  
  ## results
  name = "WSS: Weighted Sum Statistic"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(wss.stat = wss.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

