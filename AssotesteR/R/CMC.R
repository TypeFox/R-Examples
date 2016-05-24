#'CMC: Combined Multivariate and Collapsing Method
#'
#'The CMC method is a pooling approach proposed by Li and Leal (2008) that uses
#'allele frequencies to determine the partition of the variants into groups.
#'After the rare variants are selected, they are collapsed into an indicator
#'variable, and then a multivariate test such as Hotelling's T2 test is applied
#'to the collection formed by the common variants and the collapsed
#'super-variant.
#'
#'Those variants with minor allele frequency below the specified \code{maf}
#'threshold are collapsed into a single super variant \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param maf numeric value indicating the minor allele frequency threshold for
#'rare variants (\code{maf=0.05} by default)
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem cmc.stat cmc statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, rare variants, maf threshold, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{WSS}}, \code{\link{CMAT}}, \code{\link{TTEST}}
#'@references Li B, Leal SM (2008) Methods for Detecting Associations with Rare
#'Variants for Common Diseases: Application to Analysis of Sequence Data.
#'\emph{The American Journal of Human Genetics}, \bold{83}: 311-321
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
#'  phenotype = c(rep(1,cases), rep(0,controls))
#'
#'  # genotype matrix with 10 variants (random data)
#'  set.seed(1234)
#'  genotype = matrix(rbinom(total*10, 2, 0.051), nrow=total, ncol=10)
#'
#'  # apply CMC with maf=0.05 and 500 permutations
#'  mycmc = CMC(phenotype, genotype, maf=0.05, perm=500)
#'  mycmc
#'  }
#'
CMC <-
function(y, X, maf=0.05, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm

  ## number of individuals N
  N = nrow(X)
  ## get minor allele frequencies
  MAF = colMeans(X, na.rm=TRUE) / 2   
  ## how many variants < maf
  rare.maf = MAF < maf
  rare = sum(rare.maf)
  ## collapsing
  if (rare <= 1) 
  {   
    # if rare variants <= 1, then NO collapse is needed
    X.new = X
  } else {
    # collapsing rare variants into one column
    X.collaps = rowSums(X[,rare.maf], na.rm=TRUE)
    X.collaps[X.collaps != 0] = 1
    # joining collapsed to common variants
    X.new = cbind(X[,!rare.maf], X.collaps)	   
  }
  ## change values to -1, 0, 1
  X.new = X.new - 1
  ## number of new variants
  M = ncol(X.new)
  ## Hotellings T2 statistic
  cmc.stat = my_cmc_method(y, X.new)
  
  ## Asymptotic p-values
  # under the null hypothesis T2 follows an F distribution 
  f.stat = cmc.stat * (N-M-1)/(M*(N-2))
  df1 = M          # degrees of freedom  
  df2 = N - M - 1  # degrees of freedom  
  asym.pval = 1 - pf(f.stat, df1, df2)
  
  ## under the alternative hyposthesis T2 follows a chi-square distr
  # pval = 1 - pchisq(cmc.stat, df=M)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_cmc_method(y[perm.sample], X.new) 
    }
    # p-value 
    perm.pval = sum(x.perm > cmc.stat) / perm
  }
  
  ## results
  name = "CMC: Combined Multivariate and Collapsing Method"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare,  maf, perm)
  arg.spec = as.character(arg.spec)
  names(arg.spec) = c("cases", "controls", "variants", "rarevar", "maf", "perm")	
  res = list(cmc.stat = cmc.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

