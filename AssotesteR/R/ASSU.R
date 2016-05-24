#'ASSU: Adaptive Sum of Squared Score U Statistic
#'
#'The adaptive SSU test has been proposed by Han and Pan (2010) in an attempt
#'to overcome some of the drawbacks of the SSU test, by extending the idea of
#'the adaptive Neyman's test (Fan, 1996). The approach behind the adaptive test
#'is to use the U-statistics of the score test (from logistic regression
#'models) in order to construct a statistic with the first components of the
#'score vector U.
#'
#'\code{ASSU} gives the normal (unordered) test. \cr \code{ASSU.Ord} gives the
#'ordered (decreasing) test. \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@aliases ASSU ASSU.Ord
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem assu.stat assu statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SSU}}
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
#'  # apply ASSU with 500 permutations
#'  myassu = ASSU(phenotype, genotype, perm=500)
#'  myassu
#'
#'  # apply ASSU.Ord with 500 permutations
#'  myassuord = ASSU.Ord(phenotype, genotype, perm=500)
#'  myassuord
#'  }
#'
ASSU <-
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
  stat.assu = my_assu_method(U, V)
  assu.stat1 = stat.assu[1]  # stat normal
  p1.assu = stat.assu[2]     # pval normal
  
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
      perm.assu = my_assu_method(U.perm, V)
      x.perm[i] = perm.assu[2]
    }
    # p-value 
    perm.pval = sum(x.perm < p1.assu) / perm	 # normal
  }
  
  ## results
  name = "ASSU: Adaptive Sum of Squared Score"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(assu.stat = assu.stat1, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

