#'RWAS: Rare-Variant Weighted Aggregate Statistic
#'
#'The RWAS method has been proposed by Sul et al (2011) as a pooling method
#'that groups variants and computes a weighted sum of differences between case
#'and control mutation counts where weights are estimated from data. Under the
#'null hypothesis the RWAS statistic has an asymptotic standard normal
#'distribution, but a permutation procedure can also be applied to assess
#'statistical significance
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param maf numeric value indicating the minor allele frequency threshold for
#'rare variants (\code{ma f=0.05} by default)
#'@param perm positive integer indicating the number of permutations
#'(\code{NULL} by default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem rwas.stat rwas statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem perm.pval permuted p-value, only when \code{perm} is used
#'@returnItem args descriptive information with number of controls, cases,
#'variants, rare variants, maf and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{CMC}}
#'@references Sul JH, Han B, He D, Eskin E (2011) An Optimal Weighted
#'Aggregated Association Test for Identification of Rare Variants Involved in
#'Common Diseases. \emph{Genetics}, \bold{188}: 181-188
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
#'  set.seed(1234)  
#'  genotype = matrix(rbinom(total*10, 2, 0.051), nrow=total, ncol=10)
#'
#'  # apply RWAS with maf=0.05 and 500 permutations
#'  myrwas = RWAS(phenotype, genotype, maf=0.05, perm=500)
#'  myrwas
#'  }
#'
RWAS <-
function(y, X, maf=0.05, perm=NULL)
{
  ## checking argumetns
  if (!is.vector(y) || mode(y) != "numeric")
    stop("argument 'y' must be a numeric vector")
  if (any(is.na(y))) 
    stop("Sorry =(   No missing data allowed in argument 'y' ")	
  if (!all(y %in% c(0, 1)))
    stop("Sorry =(   argument 'y' must contain only 0 and 1")
  if(!is.matrix(X) & !is.data.frame(X))
    stop("argument 'X' must be a matrix or data.frame")
  if (nrow(X) != length(y)) 
    stop("'X' and 'y' have different lengths")
  if (!is.matrix(X)) X = as.matrix(X)
  if (mode(maf)!= "numeric" || length(maf) != 1 || maf<=0 || maf>1)
    stop("argument 'maf' incorreclty defined; must be a value between 0 and 1")
  
  #    if (!is.null(weights))
  #	{ 
  #	    if (mode(weights) != "numeric" || !all(weights >= 0)) 
  #		    stop("argument 'weights' must contain non-negative numbers")
  #       if (length(weights) != ncol(X)) 
  #	        stop("length of 'weights' differs from number of columns in 'X'")
  #    } else {
  #	    weights = rep(1, ncol(X))
  #	}
  
  if (!is.null(perm))
  {
    if (mode(perm) != "numeric" || length(perm) != 1
        || perm < 0 || (perm %% 1) !=0) 
    {
      warning("Argument 'perm' incorrectly defined. Value perm=100 is used")
      perm = 100
    }
  } else perm=0
  
  weights = rep(1, ncol(X))
  ## get minor allele frequencies
  MAFs = colSums(X, na.rm=TRUE) / (2*nrow(X))
  ## are there any rare variants?
  rare = sum(MAFs < maf) 
  if (rare == 0)
    stop(paste("\n", "Oops: No rare variants below maf=", 
               maf, " were detected. Try a larger maf", sep=""))
  ## get only rare variants
  X.new = X[ , MAFs < maf]	
  weights = weights[MAFs < maf]	
  
  ## running rwas
  rwas.stat = my_rwas_method(y, X.new, weights)	
  asym.pval = 1 - pnorm(rwas.stat)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)	     
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_rwas_method(y[perm.sample], X.new, weights) 
    }
    # p-value 
    perm.pval = sum(x.perm > rwas.stat) / perm
  }
  
  ## results
  name = "RWAS: Rare-Variant Weighted Aggregate Statistic"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare, maf, perm)
  arg.spec = as.character(arg.spec)
  names(arg.spec) = c("cases", "controls", "variants", "rarevar", "maf", "n.perms")	
  res = list(rwas.stat = rwas.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}
