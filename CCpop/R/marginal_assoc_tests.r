# Tests for association (or independence) between a SNP locus and a binary 
# disease phenotype.

marginal.assoc.test.unconstrained.chisq = function(t0, t1) {
  # The marginal Pearson chi-squared test on the 2x3 (or less) contingency table.
  
  # NOTE: the tables must be pruned so that they do not contain columns with all 
  # zeros.
  
  n0 = sum(t0)
  n1 = sum(t1)
  nx = t0 + t1
  n = n0 + n1
  e0 = (n0 / n) * nx
  e1 = (n1 / n) * nx
  
  statistic = sum((t0 - e0) ^ 2 / e0 + (t1 - e1) ^ 2 / e1)
  p.value = pchisq(statistic, df = length(t0), lower.tail = F)
  
  return (list(pen = t1 / (t0 + t1), statistic = statistic, p.value = p.value))
}

marginal.assoc.test.unconstrained.gsq = function(t0, t1) {
  # The marginal G^2 test (GLRT conditioned on the observed marginals of the 2x3
  # contingency table)
  
  pen.full = t1 / (t0 + t1)
  pen.null = sum(t1) / sum(t1 + t0)
  
  statistic = 2 * (sum(t0 * log(1 - pen.full) + t1 * log(pen.full), na.rm = T) - 
                   sum(t0 * log(1 - pen.null) + t1 * log(pen.null), na.rm = T))
  p.value = pchisq(statistic, df = length(t0), lower.tail = F)
  
  return (list(pen = pen.full, statistic = statistic, p.value = p.value))
}

marginal.assoc.test.hwe.in.controls = function(t0, t1) {
  # Assuming HWE among the controls (justifiable if the disease is rare). This is 
  # a Wald test for marginal association assuming HWE among the controls (Chen and
  # Chatterjee 2007)
  
  t1[t1 == 0] = 1 # is this a good solution to the case where I have zero t1 cells?
  n0 = sum(t0)
  
  ef = (t0[2] + 2 * t0[3]) / (2 * n0)
  et0 = n0 * c((1 - ef) ^ 2, 2 * ef * (1 - ef), ef ^ 2)
  eor.hom = (t1[2] * et0[1]) / (et0[2] * t1[1])
  eor.het = (t1[3] * et0[1]) / (et0[3] * t1[1])
  
  beta = log(c(eor.hom, eor.het))
  
  beta.cov = matrix(c(
    1/t1[1] + 1/t1[2] + 1/(2*et0[1] + et0[2]) + 1/(2*et0[3] + et0[2]),
    1/t1[1] + 1/n0 * 1/(ef * (1 - ef)),
    1/t1[1] + 1/n0 * 1/(ef * (1 - ef)),
    1/t1[1] + 1/t1[3] + 4/(2*et0[3] + et0[2]) + 4/(2*et0[1] + et0[2])), nrow = 2)
  
  statistic = t(beta) %*% solve(beta.cov) %*% beta
  p.value = pchisq(statistic, df = 2, lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

marginal.assoc.test.pop.hwe.kpy = function(t0, t1, tp, prevalence, pen.initial = NULL, f.initial = NULL) {
  # Assuming HWE in the general population, and that disease prevalence is 
  # known. This test uses an extraneous population sample to improve estimation 
  # of population parameters such as SNP distributions.
  
  tx = t0 + t1 + tp
  
  f.null.est = min(0.5, (0.5 * tx[2] + tx[3]) / sum(tx))
  px.null.est = c((1 - f.null.est) ^ 2, 2 * f.null.est * (1 - f.null.est), f.null.est ^ 2)
  null.loglik = sum(t0 * log(1 - prevalence) + t1 * log(prevalence) + tx * log(px.null.est))
  
  if (is.null(pen.initial)) {
    pars.initial = c(rep(prevalence, 3), f.null.est)
  } else {
    pars.initial = c(pen.initial, f.initial)
  }
  
  eval_f = function(pars) {
    pen = pars[1:3]
    f = pars[4]
    
    px = c((1 - f) ^ 2, 2 * f * (1 - f), f ^ 2)
    dpx.df = c(2 * f - 2, 2 - 4 * f, 2 * f)
    
    objective = null.loglik - sum(t0 * log(1 - pen) + t1 * log(pen)) - sum(tx * log(px))
    gradient = -c(t1 / pen - t0 / (1 - pen), sum(tx / px * dpx.df))
    
    return (list(objective = objective, gradient = gradient))
  }
  
  lb = rep(.Machine$double.eps, 4)
  ub = c(rep(1, 3) - .Machine$double.neg.eps, 0.5)
  
  eval_g_eq = function(pars) {
    pen = pars[1:3]
    f = pars[4]
    
    px = c((1 - f) ^ 2, 2 * f * (1 - f), f ^ 2)
    dpx.df = c(2 * f - 2, 2 - 4 * f, 2 * f)
    
    constraints = sum(pen * px) - prevalence
    jacobian = matrix(c(px, sum(pen * dpx.df)), nrow = 1)
    
    return (list(constraints = constraints, jacobian = jacobian))
  }

  sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
               opts = list(algorithm = 'NLOPT_LD_SLSQP', ftol_abs = 1e-8, xtol_rel = 0, maxeval = 1000))
  
  pen = t0 * 0 + sol$solution[1:3]
  f = sol$solution[4]
  px = c((1 - f) ^ 2, 2 * f * (1 - f), f ^ 2)
    
  statistic = -2 * sol$objective
  p.value = pchisq(statistic, df = 2, lower.tail = F)
  
  return (list(pen = pen, statistic = statistic, p.value = p.value))
}

marginal.assoc.test.kpx.kpy = function(t0, t1, prevalence, px, pen.initial = NULL) {
  # Assuming the SNP distribution is known, and that disease prevalence is known
  
  # The GLRT for testing a marginal association given a known SNP distribution 
  # (specified as genotypic, but can come from a HWE assumption and known MAF)
  
  # NOTE: this is exactly the same algorithm as that for a the pairwise test with 
  # known pxx.
  
  m = length(px)
  df = m - 1
  
  if (is.null(pen.initial)) {
    pars.initial = rep(prevalence, df)
  } else {
    pen.initial = pen.initial[2:m]
  }
  
  fn = function(x) {
    pen = c((prevalence - sum(x * px[2:m])) / px[1], x)
    return (-sum(t0 * log(1 - pen) + t1 * log(pen)))
  }
  
  gr = function(x) {
    x0 = (prevalence - sum(x * px[2:m])) / px[1]
    g = t1[2:m] / x - t0[2:m] / (1 - x) - t1[1] / x0 / px[1] * px[2:m] + t0[1] / (1 - x0) / px[1] * px[2:m]      
    return (-g)
  }
  
  ui = rbind(diag(df), -diag(df), -px[2:m], px[2:m])
  ci = c(rep(.Machine$double.eps, df), rep(-(1 - .Machine$double.neg.eps), df), .Machine$double.eps - prevalence, prevalence - px[1] * (1 - .Machine$double.neg.eps))
  sol = constrOptim(pars.initial, fn, gr, ui, ci)
  
  pen = t0 * 0 + c((prevalence - sum(sol$par * px[2:m])) / px[1], sol$par)
  
  statistic = 2 * sum(t0 * log((1 - pen) / (1 - prevalence)) + t1 * log(pen / prevalence))
  p.value = pchisq(statistic, df = df, lower.tail = F)
  
  return (list(pen = pen, statistic = statistic, p.value = p.value))
}
