#'CMAT: Cumulative Minor Allele Test
#'
#'CMAT is a pooling method proposed by Zawistowski et al (2010). CMAT works by
#'comparing weighted minor-allele counts (for cases and controls) against the
#'weighted major-allele counts (for cases and controls). Although the CMAT test
#'statistic is based on a chi-square statistic, it does not follow a known
#'distribution and its significance has to be determined by a permutation
#'procedure.
#'
#'By default, argument \code{maf=NULL} meaning that no rare variants are
#'selected \cr
#'
#'By default, argument \code{weights=NULL} but different values for the
#'variants can be provided \cr
#'
#'Statistical significance is determined by permutation procedure \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param maf optional numeric value to specify a threshold for the minor allele
#'frequency of rare variants (\code{NULL} by default)
#'@param weights optional vector of weights for the variants (\code{NULL} by
#'default)
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem cmat.stat cmat statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, rare variants, maf, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{CMC}}, \code{\link{WSS}}
#'@references Zawistowski M, Gopalahrishnan S, Ding J, Li Y, Grimm S, Zollner S
#'(2010) \emph{The American Journal of Human Genetics}, \bold{87}: 604-617
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
#'  # apply CMAT with 500 permutations
#'  mycmat1 = CMAT(phenotype, genotype, perm=500)
#'  mycmat1
#'
#'  # apply CMAT with maf=0.05 and 500 permutations
#'  mycmat2 = CMAT(phenotype, genotype, maf=0.05, perm=500)
#'  mycmat2
#'  }
#'
CMAT <- 
function(y, X, maf=NULL, weights=NULL, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  # check maf
  if (!is.null(maf))
  {
    if (mode(maf) != "numeric" || length(maf) != 1
        || maf <= 0  || maf >= 1)
      stop("argument 'maf' must be a value between 0 and 1")
  }
  # check weights
  if (!is.null(weights))
  { 
    if (mode(weights) != "numeric" || !all(weights >= 0)) 
      stop("argument 'weights' must contain non-negative numbers")
    if (length(weights) != ncol(X)) 
      stop("length of 'weights' differs from number of columns in 'X'")
  } else {
    weights = rep(1, ncol(X))
  }
  
  # get rare variants if maf != NULL
  rare = "NULL"
  if (!is.null(maf))
  {
    ## get minor allele frequencies
    MAFs = colMeans(X, na.rm=TRUE) / 2    
    ## are there any rare variants?
    rare = sum(MAFs < maf)
    if (rare > 0) 
    {
      X.new = X[ , MAFs < maf]
      weights = weights[MAFs < maf]
      rare = sum(weights != 0)
      if (sum(weights) == 0)
        stop(paste("Oops!, with maf =", maf, ", all weights are zero"))
    }		
    if (rare == 0)
      stop(paste("Ooops, no rare variants detected below maf=", maf, sep=""))	
  } else {
    X.new = X
    maf = "NULL"
  }
  
  ## apply cmat.method
  cmat.stat = my_cmat_method(y, X.new, weights)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_cmat_method(y[perm.sample], X.new, weights) 
    }
    ## p-value 
    perm.pval = sum(x.perm > cmat.stat) / perm
  }
  
  ## results
  name = "CMAT: Cumulative Minor Allele Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, perm)
  names(arg.spec) = c("cases", "controls", "variants", "rarevar", "maf", "n.perms")	
  res = list(cmat.stat = cmat.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

