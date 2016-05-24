#'ASCORE: Adaptive Score Test
#'
#'The Adaptive Score test has been proposed by Hand and Pan (2010) in an
#'attempt to overcome some of the drawbacks of the Score test (from logistic
#'regression) by extending the idea of the adaptive Neymans test. The approach
#'behind the adaptive test is to use the the first components of the U-score
#'vector in order to construct a test statistic
#'
#'\code{ASCORE} gives the normal (unordered) test. \cr \code{ASCORE.Ord} gives
#'the ordered (decreasing) test. \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@aliases ASCORE ASCORE.Ord
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem ascore.stat ascore statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SCORE}}
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
#'  # apply ASCORE with 500 permutations
#'  myascore = ASCORE(phenotype, genotype, perm=500)
#'  myascore
#'
#'  # apply ASCORE.Ord with 500 permutations
#'  myascoreord = ASCORE.Ord(phenotype, genotype, perm=500)
#'  myascoreord
#'  }
#'
ASCORE <-
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
  ## run score method
  stat.sco = my_ascore_method(U, V)
  score.stat = stat.sco[1]
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    ymean = mean(y)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      # center phenotype y
      y.perm = y[perm.sample] - ymean
      # get score vector
      U.perm = colSums(y.perm * X, na.rm=TRUE)
      perm.sco = my_ascore_method(U.perm, V)
      x.perm[i] = perm.sco[1]
    }
    # p-value 
    perm.pval = sum(x.perm > score.stat) / perm	
  }
  ## results
  name = "ASCORE: Adaptive Score Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(ascore.stat = score.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

