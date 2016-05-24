#'CALPHA: C-alpha Score Test
#'
#'The C-alpha score-test of Neyman and Scott (1966) has been proposed to be
#'used in association astudies by Neale et al (2011) in order to test the
#'observed distributino of rare variants in cases versis controls. Under the
#'null hypothesis of no association between the variants and the phenotype,
#'C-alpha assumes that the distribution of counts (copies of an observed
#'variant) should follow a binomial distribution. The C-alpha test statistic
#'contrasts the variance of each observed count with the expected variance,
#'assuming the binomial distribution. Under the null hypothesis, the test
#'statistic follows a standard normal distribution
#'
#'There is no imputation for missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y Numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X Numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param perm Positive integer indicating the number of permutations
#'(\code{NULL} by default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem calpha.stat c-alpha statistic
#'@returnItem asym.pval asymptotic p-value
#'@returnItem perm.pval permuted p-value; only when \code{perm} is used
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@references Neyman J, Scott E (1966) On the use of c-alpha optimal tests of
#'composite hypothesis. \emph{Bulletin of the International Statistical
#'Institute}, \bold{41}: 477-497 \cr
#'
#'Neale BM, Rivas MA, Voight BF, Altshuler D, Devlin B, Orho-Melander M,
#'Kathiresan S, Purcell SM, Roeder K, Daly MJ (2011) Testing for an unusual
#'distribution of rare variants. \emph{PLoS Genetics}, \bold{7(3)}: e1001322
#'\cr
#'
#'Basu S, Pan W (2011) Comparison of Statistical Tests for Disease Association
#'With Rare Variants. \emph{Genetic Epidemiology}, \bold{35(7)}: 606-619
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
#'  set.seed(123)
#'  genotype = matrix(rbinom(total*10, 2, 0.10), nrow=total, ncol=10)
#'
#'  # apply CALPHA with 500 permutations
#'  mycalpha = CALPHA(phenotype, genotype, perm=500)
#'
#'  # this is what we get
#'  mycalpha
#'  }
#'
CALPHA <-
function(y, X, perm=NULL)
{  
  ## checking arguments
  if (!is.null(perm))
  {
    if (mode(perm) != "numeric" || length(perm) != 1
        || perm < 0 || (perm %% 1) !=0) 
    {
      warning("Argument 'perm' incorrectly defined. Value perm=100 is used")
      perm = 100
    }
  } else perm=0
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  
  # how many cases
  nA = sum(y)
  # how many controls
  nU = sum(y == 0)
  # proportion of cases
  p0 = nA / (nA + nU)
  # how many variants
  m = ncol(X)
  # copies of the i-th variant type
  n = apply(X, 2, function(x) sum(x>0, na.rm=TRUE))
  # copies of the i-th variant type in the cases
  g = apply(X[y==1,], 2, function(x) sum(x>0, na.rm=TRUE))
  
  # Test statistic 
  calpha.stat = my_calpha_method(y, X)
  # Variance of Talpha
  Valpha = 0
  for (i in 1:m) {
    for (u in 0:n[i]) {
      Valpha = Valpha + (((u - n[i]*p0)^2 - n[i]*p0*(1-p0))^2)*dbinom(u, n[i], p0)
    }
  }
  names(Valpha) = NULL
  # Z score
  Zscore = calpha.stat / sqrt(Valpha)
  
  # asymptotic p-vaue
  if (Valpha==0) asym.pval=1 else
    asym.pval = 1 - pchisq(calpha.stat^2 / Valpha, df=1)
  
  # permutations
  perm.pval = NA
  if (perm != 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_calpha_method(y[perm.sample], X) 
    }
    # p-value 
    perm.pval = sum(x.perm^2 > calpha.stat^2) / perm
  }	  
  
  ## results
  name = "CALPHA: c-alpha Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")
  res = list(calpha.stat = calpha.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"	
  return(res)
}

