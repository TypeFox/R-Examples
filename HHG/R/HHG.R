# Main R file for the Heller-Heller-Gorfine test package


# The general test of independence (with or without handling of ties)
hhg.test = function(Dx, Dy, ties = T, w.sum = 0, w.max = 2, nr.perm = 10000, 
  is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, seq.alpha0 = NULL, 
  seq.beta0 = NULL, seq.eps = NULL, nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.double(Dx) || !is.double(Dy) || !is.matrix(Dx) || !is.matrix(Dy) || 
        nrow(Dx) != ncol(Dx) || nrow(Dx) != nrow(Dy) || nrow(Dy) != ncol(Dy)) {
    stop('Dx and Dy are expected to be square matrices of doubles, and must have the same number of rows/cols')
  }
  if (w.sum < 0 || w.max < 0) {
    stop('w.sum and w.max should be greater or equal to zero')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  
  if (ties) {
    test_type = .MV_IND_HHG
  } else {
    test_type = .MV_IND_HHG_NO_TIES
  }

  dummy.y = matrix(0, nrow(Dy), 1)
  extra_params = as.double(0)
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, seq.total.nr.tests, seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, dummy.y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  
  class(ret) = 'HHG.Test.Result'
  ret$stat.type = 'hhg.test'
  ret$n = nrow(Dx)
  return (ret)
}

# The K-sample test (with y_i in 0:(K-1))
hhg.test.k.sample = function(Dx, y, w.sum = 0, w.max = 2, nr.perm = 10000, 
  is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, seq.alpha0 = NULL, 
  seq.beta0 = NULL, seq.eps = NULL, nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.double(Dx) || !is.matrix(Dx) || 
        nrow(Dx) != ncol(Dx) || nrow(Dx) != length(y)) {
    stop('Dx is expected to be a square matrix of doubles, and must have the same number of rows/cols as the vector y')
  }
  if (w.sum < 0 || w.max < 0) {
    stop('w.sum and w.max should be greater or equal to zero')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (any(y != round(y))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (max(y) > 10) {
    warning('the K-sample test is appropriate for small values of K, consider using the general test, implemented by hhg.test')
  }
  
  test_type = .MV_KS_HHG
  
  dummy.Dy = 0
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  extra_params = as.double(0)
  is_sequential = as.integer(is.sequential)

  wald = .configure.wald.sequential(is.sequential, seq.total.nr.tests, seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
    
  res = .Call('HHG_R_C', test_type, Dx, dummy.Dy, y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  
  class(ret) = 'HHG.Test.Result'
  ret$stat.type = 'hhg.test.k.sample'
  ret$n = nrow(Dx)
  ret$y = table(y)
  return (ret)
}


# The 2-sample test (with y_i in {0, 1})
hhg.test.2.sample = function(Dx, y, w.sum = 0, w.max = 2, nr.perm = 10000, 
  is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, seq.alpha0 = NULL, 
  seq.beta0 = NULL, seq.eps = NULL, nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (!is.double(Dx) || !is.matrix(Dx) || nrow(Dx) != ncol(Dx) || nrow(Dx) != length(y)) {
    stop('Dx is expected to be a square matrix of doubles, and must have the same number of rows/cols as the vector y')
  }
  if (w.sum < 0 || w.max < 0) {
    stop('w.sum and w.max should be greater or equal to zero')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (!all(y %in% c(0, 1))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1}')
  }
  
  test_type = .MV_TS_HHG
  
  dummy.Dy = 0
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  extra_params = as.double(0)
  is_sequential = as.integer(is.sequential)

  wald = .configure.wald.sequential(is.sequential, seq.total.nr.tests, seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
    
  res = .Call('HHG_R_C', test_type, Dx, dummy.Dy, y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  
  class(ret) = 'HHG.Test.Result'
  ret$stat.type = 'hhg.test.2.sample'
  ret$n = nrow(Dx)
  ret$y = table(y)
  
  return (ret)
}



#function for printing the univariate object
print.HHG.Test.Result = function(x,...){
  if(x$stat.type == 'hhg.test'){
    cat(paste0('Results for the HHG test of independence:\n'))
    cat(paste0('Number of observations: ',x$n,'\n'))
  }else if(x$stat.type == 'hhg.test.2.sample'){
    cat(paste0('Results for the HHG test for equality of distributions, with two groups:\n'))
    cat(paste0('Number of observations: ',x$n,'\n'))
    cat(paste0('Group sizes: \n'))
    cat(paste(x$y))
    cat('\n')
  }else if(x$stat.type == 'hhg.test.k.sample'){
    cat(paste0('Results for the HHG test for equality of distributions, with K>=2 groups:\n'))
    cat(paste0('Number of observations: ',x$n,'\n'))
    cat(paste0('Group sizes: \n'))
    cat(paste(x$y))
    cat('\n')
  }
  cat('\n')
  cat(paste0('Sum of chi-squared scores of all 2X2 tables: ',format(x$sum.chisq,nsmall = 3),'\n'))
  if('perm.pval.hhg.sc' %in% names(x)){
    cat(paste0('P-value for the sum of chi-squared scores statistic : ',format(x$perm.pval.hhg.sc,digits = 3),'\n'))
  }
  cat('\n')
  
  cat(paste0('Sum of likelihood ratio scores of all 2X2 tables: ',format(x$sum.lr,nsmall = 3),'\n'))
  if('perm.pval.hhg.sl' %in% names(x)){
    cat(paste0('P-value for the sum of likelihood ratio scores statistic : ',format(x$perm.pval.hhg.sl,digits = 3),'\n'))
  }
  cat('\n')
  
  cat(paste0('Maximum over chi-squared scores of all 2X2 tables: ',format(x$max.chisq,nsmall = 3),'\n'))
  if('perm.pval.hhg.mc' %in% names(x)){
    cat(paste0('P-value for the max of chi-squared scores statistic : ',format(x$perm.pval.hhg.mc,digits = 3),'\n'))
  }
  cat('\n')
  
  cat(paste0('Maximum over likelihood ratio scores of all 2X2 tables: ',format(x$max.lr,nsmall = 3),'\n'))
  if('perm.pval.hhg.ml' %in% names(x)){
    cat(paste0('P-value for the max of likelihood ratio scores statistic : ',format(x$perm.pval.hhg.ml,digits = 3),'\n'))
  }
  cat('\n')
}



