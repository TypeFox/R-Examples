#'RBT: Replication Based Test
#'
#'The Replication-Based Test has been proposed by Ionita-Laza et al (2011). RBT
#'is based on a weighted-sum statistic that is similar to a pooled association
#'test but purposefully designed to deal with possibly different association
#'directions. The significance of the statistic has to be calculated by
#'permutations.
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
#'@returnItem rbt.stat rbt statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{WSS}}, \code{\link{CMC}}
#'@references Ionita-Laza I, Buxbaum JD, Laird NM, Lange C (2011) A New Testing
#'Strategy to Identify Rare Variants with Either risk or Protective Effects on
#'Disease. \emph{PLoS Genetics}, \bold{7(2)}: e1001289
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
#'  # apply RBT with 500 permutations
#'  myrbt = RBT(phenotype, genotype, perm=500)
#'  myrbt
#'  }
#'
RBT <-
function(y, X, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  	
  ## apply RBT method
  rbt.stat = my_rbt_method(y, X)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_rbt_method(y[perm.sample], X) 
    }
    # p-value 
    perm.pval = sum(x.perm > rbt.stat) / perm
  }
  
  ## results
  name = "RBT: Replication Based Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(rbt.stat = rbt.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

