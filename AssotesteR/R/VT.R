#'VT: Variable Threshold
#'
#'The Variable Threshold (VT) test has been proposed by Price et al (2010)
#'based on the assumption that the minor allele frequencies of the causal rare
#'variants may be different from those nonfunctional rare variants, which, if
#'true, can be utilized to improve the power of the corresponding pooled
#'association tests.  The idea behind this approach is that there exists some
#'(unknown) threshold T for which variants with a minor allele frequency (MAF)
#'below T are more likely to be functional than are variants with an MAF above
#'T. VT works by finding the maximum z-score across all possible values for the
#'threshold T.
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param maf numeric value indicating the minor allele frequency threshold for
#'rare variants (must be a positive number between 0 and 1, \code{maf=0.05} by
#'default)
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem vt.stat vt statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{WSS}}
#'@references Price AL, Kryukov GV, de Bakker PIW, Purcell SM, Staples J, Wei
#'LJ, Sunyaev SR (2010) Pooled Association Tests for Rare Variants in
#'Exon-Sequencing Studies. \emph{The American Journal of Human Genetics},
#'\bold{86}: 832-838
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
#'  #  phenotype vector
#'  phenotype = c(rep(1, cases), rep(0, controls))
#'
#'  # genotype matrix with 10 variants (random data)
#'  set.seed(1234)  
#'  genotype = matrix(rbinom(total*10, 2, 0.051), nrow=total, ncol=10)
#'
#'  # apply VT with maf=0.05 and 500 permutations
#'  myvt = VT(phenotype, genotype, maf=0.05, perm=500)
#'  myvt
#'  }
#'
VT <-
function(y, X, maf=0.05, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  if (mode(maf) != "numeric" || length(maf) != 1
      || maf <= 0  || maf >= 1)
    stop("argument 'maf' must be a value between 0 and 1")
  
  ## running vt method
  mafs = (1 + colSums(X, na.rm=TRUE)) / (2 + 2*nrow(X))
  h.maf = sort(unique(mafs))
  vt.stat = my_vt_method(y, X, mafs, h.maf)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_vt_method(y[perm.sample], X, mafs, h.maf)
    }
    ## p-value
    perm.pval = sum(x.perm >= vt.stat) / perm
  }
  
  ## results
  name = "VT: Variable Threshold"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), maf, perm)
  names(arg.spec) = c("cases", "controls", "variants", "maf", "n.perms")	
  res = list(vt.stat = vt.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

