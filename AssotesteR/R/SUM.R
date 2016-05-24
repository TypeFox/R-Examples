#'SUM: Sum Test
#'
#'The SUM test has been proposed by Pan (2009) based on a modification of the
#'Score. The idea behind the Sum Test is to test on only one parameter under
#'the assumption of a common association strength between each of multiple
#'genetic variants (e.g. SNPs) and the trait under analysis. The Sum test
#'focuses on a scalar function of the multiple parameters with a resulting
#'degree of freedom DF=1
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
#'@returnItem sum.stat sum statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SCORE}}, \code{\link{SSU}}, \code{\link{SSUW}},
#'\code{\link{WST}}
#'@references Pan W (2009) Asymptotic tests of association with multiple SNPs
#'in linkage disequilibrium. \emph{Genetic Epidemiology}, \bold{33}: 497-507
#'\cr
#'
#'Pan W, Han F, Shen X (2010) Test Selection with Application to Detecting
#'Association with Multiple SNPs. \emph{Human Heredity}, \bold{69}: 120-130
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
#'  # apply SUM with 500 permutations
#'  mysum = SUM(phenotype, genotype, perm=500)
#'  mysum
#'  }
#'
SUM <-
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
  ## run sum method
  stat.sum = my_sum_method(U, V)
  sum.stat = stat.sum[1]
  asym.pval = stat.sum[2]
  
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
      perm.sco = my_sum_method(U.perm, V)
      x.perm[i] = perm.sco[1]
    }
    # p-value 
    perm.pval = sum(x.perm > sum.stat) / perm	
  }
  
  ## results
  name = "SUM: Sum Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(sum.stat = sum.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

