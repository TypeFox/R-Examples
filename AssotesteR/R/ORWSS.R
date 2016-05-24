#'ORWSS: Odds Ratio Weighted Sum Statistic
#'
#'The ORWSS method has been proposed by Feng et al (2011) and it is based on a
#'weighted sum statistic like the WSS method of Madsen and Browning (2009).
#'ORWSS uses the logarithm of the odds ratio of a genetic variant as the weight
#'for that variant, rather than the variance estimated in controls.
#'
#'When \code{c.param=NULL}, the weights of the sum statistic are simply the
#'logarithm of the amended Odds Ratio of each variant (as in Dai et al 2012).
#'Alternative values like \code{c.param=1.64} or \code{c.param=1.28} are
#'suggested in Feng et al (2011). \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param c.param optional value to specify the \code{c} parameter. See
#'reference Feng et al, 2011
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem orwss.stat orwss statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{WSS}}
#'@references Feng T, Elston RC, Zhu X (2011) Detecting Rare and Common
#'Variants for Complex Traits: Sibpair and Odds Ratio Weighted Sum Statistics
#'(SPWSS, ORWSS). \emph{Genetic Epidemiology}, \bold{35}: 398-409 \cr
#'
#'Dai Y, Jiang R, Dong J (2012) Weighted selective collapsing strategy for
#'detecting rare and common variants in genetic association study. \emph{BMC
#'Genetics}, \bold{13}:7
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
#'  # apply ORWSS with c.param=NULL and 500 permutations
#'  myorwss1 = ORWSS(phenotype, genotype, c.param=NULL, perm=100)
#'  myorwss1
#'
#'  # apply ORWSS with c.param=1.64 (see Feng et al 2011)
#'  myorwss2 = ORWSS(phenotype, genotype, c.param=1.64, perm=100)
#'  myorwss2
#'  }
#'
ORWSS <- 
function(y, X, c.param=NULL, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  if (!is.null(c.param))
  {
    if (mode(c.param)!= "numeric" || length(c.param) != 1
        || c.param<0 || c.param>=10)
      stop("argument 'c.param' incorreclty defined; must be a non-negative value")
  }
  
  Xnew = X
  Xnew[Xnew!=0] = 1
  ## running orwss method
  orwss.stat = my_orwss_method(y, Xnew, c.param)	
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_orwss_method(y[perm.sample], Xnew, c.param) 
    }
    ## p-value 
    perm.pval = sum(x.perm > orwss.stat) / perm
  }
  
  ## results
  name = "ORWSS: Odds Ratio Weighted Sum Statistic"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(orwss.stat = orwss.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}





