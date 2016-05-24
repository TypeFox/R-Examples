#'ASSUW: Adaptive Sum of Squared Score U Statistic
#'
#'The adaptive Weighted Score test has been proposed by Han and Pan (2010) in
#'an attempt to overcome some of the drawbacks of the SSUW test, by extending
#'the idea of the adaptive Neyman's test (Fan, 1996). The approach behind the
#'adaptive test is to use the U-statistics of the score test (from logistic
#'regression models) in order to construct a statistic with the first
#'components of the score vector U.
#'
#'\code{ASSUW} gives the normal (unordered) test. \cr \code{ASSUW.Ord} gives
#'the ordered (decreasing) test. \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@aliases ASSUW ASSUW.Ord
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem assuw.stat assuw statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SSUW}}
#'@references Han F, Pan W (2010) A Data-Adaptive Sum Test for Disease
#'Association with Multiple Common or Rare Variants. \emph{Human Heredity},
#'\bold{70}: 42-54 \cr
#'
#'Pan W, Shen X (2011) Adaptive Tests for Association of Rare Variants.
#'\emph{Genetic Epidemiology}, \bold{35}: 381-388
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
#'  # apply ASSUW with 500 permutations
#'  myassuw = ASSUW(phenotype, genotype, perm=500)
#'  myassuw
#'
#'  # apply ASSUW.Ord with 500 permutations
#'  myassuword = ASSUW.Ord(phenotype, genotype, perm=500)
#'  myassuword
#'  }
#'
ASSUW <-
function(y, X, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  
  ## get U and V
  getuv = my_getUV(y, X)
  U = getuv$U
  V = getuv$V
  ## run assuw method
  stat.assuw = my_assuw_method(U, V)
  assuw.stat1 = stat.assuw[1]  # stat normal
  p1.assuw = stat.assuw[2]     # pval normal
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    p1.perm = rep(0, perm)
    ymean = mean(y)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      # center phenotype y
      y.perm = y[perm.sample] - ymean
      # get score vector
      U.perm = colSums(y.perm * X, na.rm=TRUE)
      perm.assuw = my_assuw_method(U.perm, V)
      p1.perm[i] = perm.assuw[2]
    }
    # p-value 
    perm.pval = sum(p1.perm < p1.assuw) / perm  # normal
  }
  
  ## results
  name = "ASSUW: Adaptive Weighted Sum of Squared Score"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(assuw.stat = assuw.stat1, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

