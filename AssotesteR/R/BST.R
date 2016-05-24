#'BST: Bayesian Score Test
#'
#'BST is based on the test statistic of Goeman et al (2005) following a general
#'empirical Bayes method to test on a large number of parameters in a logistic
#'regression model. BST is closely related to the usual Score test, although it
#'assumes an empirical Bayesian model with an independent prior on the genetic
#'variant effects. The null distribution of the BST statistic is unknown and
#'has to be estimated by permutation.
#'
#'The BST statistic does not offer an asymptotic p-value. Permutations are
#'required \cr
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
#'@returnItem bst.stat bst statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, and permutations
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{SCORE}}
#'@references Goeman JJ, van de Geer SA, van Houwelingen HC (2006) Testing
#'against a high dimensional alternative. \emph{Journal of the Royal
#'Statistical Society}, \bold{68}: 477-493 \cr
#'
#'Chapman J, Whittaker J (2008) Analysis of Multiple SNPs in a Candidate Gene
#'or Region. \emph{Genetic Epidemiology}, \bold{32}: 560-566 \cr
#'
#'Pan W (2009) Asymptotic Tests of Association with Multiple SNPs in Linkage
#'Disequilibrium. \emph{Genetic Epidemiology}, \bold{33}: 497-507
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
#'  genotype = matrix(rbinom(total*10, 2, 0.05), nrow=total, ncol=10)
#'  
#'  # apply BST with 500 permutations
#'  mybst = BST(phenotype, genotype, perm=500)  
#'  mybst
#'  }
#'
BST <-
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
  bst.stat = my_bst_method(U, V)
  
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
      perm.bst = my_bst_method(U.perm, V)
      x.perm[i] = perm.bst
    }
    # p-value 
    perm.pval = sum(x.perm > bst.stat) / perm	
  }
  
  ## results
  name = "BST: Bayesian Score Test"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")
  res = list(bst.stat = bst.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"	
  return(res)
}

