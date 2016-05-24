#'SCORE: Score Test (from Logistic Regression)
#'
#'The Score test is one of the statistical tests used for logistic regression
#'models, which is one of the standard approaches used in genetic association
#'studies. Under the null hypothesis that there is no associated variants
#'within the examined region, the test statistic has an asymptotic chi-square
#'distribution. In addition, a permutation procedure can be applied to assess
#'the significance of the test.
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
#'@returnItem score.stat score statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SSU}}, \code{\link{SSUW}}, \code{\link{SUM}}
#'@references Clayton D, Chapman J, Cooper J (2004) Use of unphased multilocus
#'genotype data in indirect association studies. \emph{Genetic Epidemiology},
#'\bold{27}: 415-428 \cr
#'
#'Chapman J, Whittaker J (2008) Analysis of Multiple SNPs in a Candidate Gene
#'or Region. \emph{Genetic Epidemioloy}, \bold{32}: 560-566
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
#'  # phenotype (first column of genodata)
#'  phenotype = c(rep(1, cases), rep(0, controls))
#'
#'  # genotype matrix with 10 variants (random data)
#'  set.seed(123)
#'  genotype = matrix(rbinom(total*10, 2, 0.05), nrow=total, ncol=10)
#'
#'  # apply SCORE with 500 permutations
#'  myscore = SCORE(phenotype, genotype, perm=500)
#'  myscore
#'  }
#'
SCORE <-
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
  stat.sco = my_score_method(U, V)
  score.stat = stat.sco[1]
  asym.pval = stat.sco[2]
  
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
      perm.sco = my_score_method(U.perm, V)
      x.perm[i] = perm.sco[1]
    }
    # perm p-value 
    perm.pval = sum(x.perm > score.stat) / perm	
  }
  
  ## results
  name = "SCORE: Score Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(score.stat = score.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

