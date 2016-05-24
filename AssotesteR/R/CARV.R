#'CARV: Comprehrensive Approach to Analyzing Rare Variants
#'
#'The CARV method has been proposed by Hoffmann et al (2010) as an approach
#'that determines an optimal grouping of rare variants while avoiding
#'assumptions required by other methods for grouping such variants. The idea
#'behind CARV is to try multiple models for rare variants, since prior
#'information is generally not very accurate. Statistical significance is
#'obtained by permutation.
#'
#'The argument \code{waf} is used to specify weights of the variants in order
#'to incorporate allele frequency information. When \code{waf=FALSE}, all
#'variants have a constant unit weight. When \code{waf=TRUE}, the weights are
#'calculated as in the function \code{\link{WSS}}, that is, weights are the
#'inverse variance of allele frequency in controls. \cr
#'
#'The argument \code{signs} is used to specify the direction of the variant
#'effect (deleterious or protective). When \code{signs=FALSE}, all variants
#'have a positive sign indicating a likely deleterious effect. When
#'\code{signs=TRUE}, all rare alleles that are more prevalent in controls than
#'cases will have an associated sk=-1 (see paper in the reference); conversely,
#'all the rare alleles that are more prevalent in cases than controls will have
#'an associated sk=1. \cr
#'
#'The argument \code{approach} is used to specify whether the alle belongs in
#'the model for variable selection. When \code{approach="hard"}, only variants
#'below the predefined \code{maf} are included in the analysis. When
#'\code{approach="variable"}, a variable threshold approach is used as in the
#'\code{\link{VT}} method: all possible minor allele frequencies are
#'considered, selecting the maximum statistic among all of them. When
#'\code{approach="stepup"}, the step-up strategy described in Hoffman et al
#'(2010) is applied. \cr
#'
#'There is no imputation for the missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param waf logical value to indicate whether weights for the variants should
#'be calculated as in the \code{WSS} method (default \code{waf=FALSE}). See
#'deatils below
#'@param signs logical value to indicate whether signs for the variants should
#'be calculated based on allele prevalence (\code{signs=FALSE} by default). See
#'details below
#'@param approach character string to indicate the type of approach to be used
#'for variable selection; i.e. whether each variant belongs in the model for
#'variable selection (\code{approach="hard"} by default). Possible options are
#'\code{"hard"}, \code{"variable"}, and \code{"stepup"}. See details below
#'@param maf numeric value between 0 and 1 to indicate the threshold of minor
#'allele frequency for rare variants. Only used when \code{approach="hard"}
#'@param perm positive integer indicating the number of permutations (100 by
#'default)
#'@return An object of class \code{"assoctest"}, basically a list with the
#'following elements:
#'@returnItem carv.stat carv statistic
#'@returnItem perm.pval permuted p-value
#'@returnItem args descriptive information with number of controls, cases,
#'variants, permutations, waf, signs, approach, and maf
#'@returnItem name name of the statistic
#'@author Gaston Sanchez
#'@seealso \code{\link{RARECOVER}}
#'@references Hoffmann TJ, Marini NJ, Witte JS (2010) Comprehensive Approach to
#'Analyzing Rare Genetic Variants. \emph{PLoS One}, \bold{5(11)}: e13584
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
#'  # define genotype matrix with 10 variants (random data)
#'  set.seed(1234)
#'  genotype = matrix(rbinom(total*10, 2, 0.051), nrow=total, ncol=10)
#'
#'  # apply CARV with "hard" approach and maf=0.05
#'  mycarv1 = CARV(phenotype, genotype, waf=FALSE, signs=FALSE, 
#'     approach="hard", maf=0.05, perm=500)
#'  mycarv1
#'
#'  # apply CARV with "variable" approach and waf=TRUE
#'  mycarv2 = CARV(phenotype, genotype, waf=TRUE, signs=FALSE,
#'     approach="variable", perm=500)
#'  mycarv2
#'
#'  # apply CARV with "stepup" approach, waf=TRUE, and signs=TRUE
#'  mycarv3 = CARV(phenotype, genotype, waf=TRUE, signs=TRUE,
#'     approach="stepup", perm=500)
#'  mycarv3
#'
#'  }
#'
CARV <- 
function(y, X, waf=FALSE, signs=FALSE, approach="hard", maf=0.05, perm=100)
{

  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm  
  if (!is.logical(waf))
    stop("argument 'waf' must be TRUE or FALSE")
  if (!is.logical(signs))
    stop("argument 'signs' must be TRUE or FALSE") 
  if (!(approach %in% c("stepup", "variable", "hard")))
    stop("argument 'approach' must be one of 'stepup', 'variable' or 'hard'")
  if (mode(maf) != "numeric" || length(maf) != 1 || maf <= 0 || maf > 1)
    maf = 0.05
	
  # how many variants
  M = ncol(X)
  # weights 'ak' to incorporate allele frequencies
  ak = rep(1, M)
  if (waf) ak = my_weights_wss(y, X)   # weights a la madsen & browning
  
  ## signs 'sk' of the variant effect 
  sk = rep(1, M)
  if (signs)
  {
    nAs = apply(X[y==1,], 2, function(x) sum(!is.na(x)))
    nUs = apply(X[y==0,], 2, function(x) sum(!is.na(x)))	 
    maf.A = colSums(X[y==1,], na.rm=TRUE) / (2*nAs)
    maf.U = colSums(X[y==0,], na.rm=TRUE) / (2*nUs)
    # more prevalent in cases than controls
    sk[maf.U > maf.A] = -1 
  }
  
  ## prepare data for approaches
  ymean = mean(y)
  y.cen = y - ymean
  X.cen = scale(X, scale=FALSE)
  w = ak * sk	
  
  ## hard approach
  if (approach == "hard") 
  {
    # get minor allele frequencies
    MAFs = colMeans(X, na.rm=TRUE) / 2
    # are there any rare variants?
    if (sum(MAFs < maf) == 0)
      stop(paste("\n", "Ooops: No rare variants below maf=", 
                 maf, " were detected. Try a larger maf", sep=""))
    carv.stat = my_hard_approach(y.cen, X.cen, w, MAFs, maf)
  } 
  ## variable approach
  if (approach == "variable")
  {
    # get minor allele frequencies
    MAFs = colMeans(X, na.rm=TRUE) / 2
    several.maf = sort(unique(MAFs))
    carv.stat = my_variable_approach(y.cen, X.cen, w, MAFs, several.maf)
  }
  ## step-up approach
  if (approach == "stepup") { 
    carv.stat = my_stepup_approach(y.cen, X.cen, w)
  }
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      # center phenotype y
      y.perm = y[perm.sample] - ymean
      # get score vector
      x.perm[i] = switch(approach, 
                         "hard" = my_hard_approach(y.perm, X.cen, w, MAFs, maf),
                         "variable" = my_variable_approach(y.perm, X.cen, w, MAFs, several.maf),
                         "stepup" = my_stepup_approach(y.perm, X.cen, w)
                         )
    }
    # p-value 
    perm.pval = sum(x.perm >= carv.stat) / perm	
  }
  
  ## results
  if (waf) mywaf="TRUE" else mywaf="FALSE"
  if (signs) mysigns="TRUE" else mysigns="FALSE"
  name = "CARV: Comprehensive Approach to Analyzing Rare Genetic Variants"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm, mywaf, mysigns, approach, maf)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms", 
                      "waf", "signs", "approach", "maf")	
  res = list(carv.stat = carv.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

