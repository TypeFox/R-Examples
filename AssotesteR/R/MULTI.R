#'MULTI: Multiple Tests
#'
#'Performs multiple association tests
#'
#'The available tests are: \code{"WSS", "ORWSS", "RWAS", "CMC", "CMAT",
#'"CALPHA", "RBT", "SCORE", "SUM", "SSU", "SSUW", "UMINP", "BST",
#'"WST", "RVT1", "RVT2", "VT"}
#'
#'There is no imputation for missing data. Missing values are simply
#'ignored in the computations.
#'
#'@param y Numeric vector with phenotype status: 0=controls, 1=cases. No
#'missing data allowed
#'@param X Numeric matrix or data frame with genotype data coded as 0, 1, 2.
#'Missing data is allowed
#'@param tests character vector with names of the tests to be applied
#'@param maf numeric value indicating the minor allele frequency threshold for
#'rare variants (\code{maf=0.05} by default)
#'@param perm Positive integer indicating the number of permutations
#'(\code{100} by default)
#'@param weights optional vector of weights for the variants (\code{NULL} by
#'default)
#'@param c.param Optional value to specify the \code{c} parameter
#'when applying \code{ORWSS}
#'@return A data frame with test statistics and permutated p-values
#'@author Gaston Sanchez
#'@examples
#'
#'  \dontrun{
#'   
#'  # number of cases
#'  cases = 250
#'
#'  # number of controls
#'  controls = 250
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
#'  # apply MULTI with 100 permutations
#'  mymulti1 = MULTI(phenotype, genotype, c("BST", "CMC", "RWAS"), perm=100)
#'
#'  # this is what we get
#'  mymulti1
#'
#'  # create list with the following tests
#'  test_list = c("BST", "CMC", "CMAT", "CALPHA", "ORWSS", "RWAS",
#'      "RBT", "SCORE", "SUM", "SSU", "SSUW", "UMINP", "WSS", "WST")
#'
#'  # apply MULTI with 100 permutations
#'  mymulti2 = MULTI(phenotype, genotype, test_list, perm=100)
#'
#'  # this is what we get
#'  mymulti2
#'  }
#'
MULTI <- 
function(y, X, tests, maf=0.05, perm=100, weights=NULL, c.param=NULL)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  # tests
  mytests = c("WSS", "ORWSS", "RWAS", "CMC", "CMAT", "CALPHA", "RBT", 
              "SCORE", "SUM", "SSU", "SSUW", "UMINP", "BST", "WST", "RVT1", "RVT2", "VT")
  test.match = tests %in% mytests  
  if (sum(test.match) != length(tests))
    stop(paste("\n", "Sorry =(   I can't work with these tests: ", tests[!test.match]))
  if (any(duplicated(tests)))
    tests = unique(tests)
  if ("WST" %in% tests)
  {
    if (sum(complete.cases(X)) != nrow(X))
      stop("Sorry,  no missing data allowed in 'X' when using test 'WST'")
  }  
  # maf
  if (mode(maf)!= "numeric" || length(maf) != 1 || maf<=0 || maf>1)
    stop("argument 'maf' incorreclty defined; must be a value between 0 and 1")
  # weights
  if (!is.null(weights))
  { 
    if (mode(weights) != "numeric" || !all(weights >= 0)) 
      stop("argument 'weights' must contain non-negative numbers")
    if (length(weights) != ncol(X)) 
      stop("length of 'weights' differs from number of columns in 'X'")
  } else {
    weights = rep(1, ncol(X))
  }
  # c.param
  if (!is.null(c.param))
  {
    if (mode(c.param)!= "numeric" || length(c.param) != 1
        || c.param<0 || c.param>=10)
      stop("argument 'c.param' incorreclty defined; must be a non-negative value")
  }
  
  ## are there any rare variants?
  MAFs = colMeans(X, na.rm=TRUE) / 2    
  RV = MAFs < maf
  rare = sum(RV)
  if (rare > 0) 
  {
    if (sum(weights[MAFs < maf]) == 0)
      stop(paste("Oops!, all weights are zero with maf =", maf))
  }		
  
  ## mean phenotype
  ymean = mean(y)
  
  ## is there any score-type test?
  st = tests %in% c("SCORE", "SUM", "SSU", "SSUW", "UMINP", "BST", "WST")
  if (sum(st) > 0) {
    getuv = my_getUV(y, X)
    U = getuv$U
    V = getuv$V
  }
  
  ## is CMC requested?
  if ("CMC" %in% tests)
  {
    ## collapsing
    if (rare <= 1) 
    {   # NO collapse is needed
      X.new = X
    } else {
      X.collaps = rowSums(X[,RV], na.rm=TRUE)
      X.collaps[X.collaps != 0] = 1
      X.new = cbind(X[,!RV], X.collaps)	   
    }
    ## change values to -1, 0, 1
    X.new = X.new - 1
  }
  ## is RWAS requested?
  if ("RWAS" %in% tests)
  {
    X.rwas = X[,RV]
    w.rwas = rep(1, ncol(X.rwas))
  }
  ## is RVT1 or RVT2 requested?	
  if ("RVT1" %in% tests || "RVT2" %in% tests)
  {
    X.rvt = X[, RV]	
    X.rvt[X.rvt == 2] = 1 
  }
  ## is VT requested?
  if ("VT" %in% tests)
  {
    mafs = (1 + colSums(X, na.rm=TRUE)) / (2 + 2*nrow(X))
    h.maf = sort(unique(mafs))
  }
  
  TESTS = tests
  k = length(tests)
  tests = tolower(TESTS)
  test.stats = rep(0, k)
  for (i in 1:k)
  {
    aux = switch(tests[i],
                 "calpha" = my_calpha_method(y, X),
                 "cmat" = my_cmat_method(y, X[,RV], weights[RV]),  
                 "cmc" = my_cmc_method(y, X.new),   
                 "rwas" = my_rwas_method(y, X.rwas, w.rwas),
                 "wss" = my_wss_method(y, X),
                 "rbt" = my_rbt_method(y, X),
                 "orwss" = my_orwss_method(y, X, c.param),
                 "score" = my_score_method(U, V),
                 "sum" = my_sum_method(U, V),
                 "ssu" = my_ssu_method(U, V),
                 "ssuw" = my_ssuw_method(U, V),
                 "uminp" = my_uminp_method(U, V),
                 "bst" = my_bst_method(U, V),
                 "wst" = my_wst_method(U, V),
                 "rvt1" = my_rvt1_method(y - ymean, X.rvt, ymean),
                 "rvt2" = my_rvt2_method(y - ymean, X.rvt, ymean),
                 "vt" = my_vt_method(y, X, mafs, h.maf)
                 )
    test.stats[i] = aux[1]
  }
  
  ## permutations
  N = length(y)
  perm.stats = matrix(0, perm, k)
  for (p in 1:perm)
  {
    perm.sample = sample(1:N)
    # is there any score-type test?
    if (sum(st) > 0) {
      # center phenotype y
      y.perm = y[perm.sample] - ymean
      U.perm = colSums(y.perm * X, na.rm=TRUE)
    }
    # for each test
    for (j in 1:k)
    {
      aux = switch(tests[j],
                   "calpha" = my_calpha_method(y[perm.sample], X), 
                   "cmat" = my_cmat_method(y[perm.sample], X[,RV], weights[RV]), 
                   "cmc" = my_cmc_method(y[perm.sample], X.new), 
                   "rwas" = my_rwas_method(y[perm.sample], X.rwas, w.rwas), 
                   "wss" = my_wss_method(y[perm.sample], X), 
                   "rbt" = my_rbt_method(y[perm.sample], X), 
                   "orwss" = my_orwss_method(y[perm.sample], X, c.param), 
                   "score" = my_score_method(U.perm, V), 
                   "sum" = my_sum_method(U.perm, V), 
                   "ssu" = my_ssu_method(U.perm, V), 
                   "ssuw" = my_ssuw_method(U.perm, V), 
                   "uminp" = my_uminp_method(U.perm, V),
                   "bst" = my_bst_method(U.perm, V), 
                   "wst" = my_wst_method(U.perm, V),
                   "rvt1" = my_rvt1_method(y.perm, X.rvt, ymean),
                   "rvt2" = my_rvt2_method(y.perm, X.rvt, ymean),
                   "vt" = my_vt_method(y.perm, X, mafs, h.maf)
                   )
      perm.stats[p, j] = aux[1]
    }
  }
  ## perm pvals
  perm.pvals = rep(0, k)
  for (j in 1:k)
    perm.pvals[j] = sum(perm.stats[,j] > test.stats[j]) / perm
  
  ## results
  res = data.frame(statistic=test.stats, pvalue=perm.pvals)
  rownames(res) = TESTS		
  res
}
