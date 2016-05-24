#'UMINP: Univariate minP (minimum p-value)
#'
#'UMINP is the Univariate minP that tests on each single genetic variant (e.g.
#'SNP) one-by-one and then takes the minimum of their p-values, Its null
#'distribution is based on numerical integration with respect to a multivariate
#'normal distribution.
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
#'@returnItem uminp.stat uminp statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SCORE}}, \code{\link{SUM}}
#'@references Pan W (2009) Asymptotic tests of association with multiple SNPs
#'in linkage disequilibrium. \emph{Genetic Epidemiology}, \bold{33}: 497-507
#'\cr
#'
#'Pan W, Han F, Shen X (2010) Test Selection with Application to Detecting
#'Disease Association with Multiple SNPs. \emph{Human Heredity}, \bold{69}:
#'120-130
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
#'  # apply UMINP with 500 permutations
#'  myuminp = UMINP(phenotype, genotype, perm=500)
#'  myuminp
#'  }
#'
UMINP <-
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
  stat.umin = my_uminp_method(U, V)
  uminp.stat = stat.umin[1]
  asym.pval = stat.umin[2]
  
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
      perm.umin = my_uminp_method(U.perm, V)
      x.perm[i] = perm.umin[1]
    }
    # permuted p-value 
    perm.pval = sum(x.perm > uminp.stat) / perm	
  }
  
  ## results
  name = "UMINP: U minimum P-value"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(uminp.stat = uminp.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

