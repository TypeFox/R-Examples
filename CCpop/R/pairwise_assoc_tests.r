# Tests for association (or independence) between a pair of SNP loci and a 
# binary disease phenotype.

.experiment = 0

.logit = function(p) {
  return (log(p / (1 - p)))
}
    
.logistic = function(x) {
  return (1 / (1 + exp(-x)))
}

.my.chisq.test.2d = function(tbl) {
  # Prune zero rows and columns
  idx.r = (rowSums(tbl) > 0)
  idx.c = (colSums(tbl) > 0)
  tbl = tbl[idx.r, ]
  tbl = tbl[, idx.c]
  
  # Compute the statistic and p-value
  n = sum(tbl)
  tbl.e = rowSums(tbl) %o% colSums(tbl) / n
  
  statistic = sum((tbl - tbl.e) ^ 2 / tbl.e, na.rm = T)
  p.value = pchisq(statistic, df = (nrow(tbl) - 1) * (ncol(tbl) - 1), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

.my.g.test.2d = function(tbl) {
  # Prune zero rows and columns
  idx.r = (rowSums(tbl) > 0)
  idx.c = (colSums(tbl) > 0)
  tbl = tbl[idx.r, ]
  tbl = tbl[, idx.c]
  
  # Compute the statistic and p-value
  n = sum(tbl)
  tbl.e = rowSums(tbl) %o% colSums(tbl) / n
  
  statistic = 2 * sum(tbl * log(tbl / tbl.e), na.rm = T)
  p.value = pchisq(statistic, df = (nrow(tbl) - 1) * (ncol(tbl) - 1), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

.my.chisq.test.3d = function(tbl) {
  # Prune zero vectors in any dimension
  idx.x = (apply(tbl, 1, sum) > 0)
  idx.y = (apply(tbl, 2, sum) > 0)
  idx.z = (apply(tbl, 3, sum) > 0)
  tbl = tbl[idx.x, , ]
  tbl = tbl[, idx.y, ]
  tbl = tbl[, , idx.z]
  
  # Compute the statistic and p-value
  n = sum(tbl)
  tbl.e = ((apply(tbl, 1, sum) %o% apply(tbl, 2, sum)) %o% apply(tbl, 3, sum)) / (n^2)
  
  statistic = sum((tbl - tbl.e) ^ 2 / tbl.e, na.rm = T)
  p.value = pchisq(statistic, df = prod(dim(tbl)) - 1 - sum(dim(tbl) - 1), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

.my.g.test.3d = function(tbl) {
  # Prune zero vectors in any dimension
  idx.x = (apply(tbl, 1, sum) > 0)
  idx.y = (apply(tbl, 2, sum) > 0)
  idx.z = (apply(tbl, 3, sum) > 0)
  tbl = tbl[idx.x, , ]
  tbl = tbl[, idx.y, ]
  tbl = tbl[, , idx.z]
  
  # Compute the statistic and p-value
  n = sum(tbl)
  tbl.e = ((apply(tbl, 1, sum) %o% apply(tbl, 2, sum)) %o% apply(tbl, 3, sum)) / (n^2)
  
  statistic = 2 * sum(tbl * log(tbl / tbl.e), na.rm = T)
  p.value = pchisq(statistic, df = prod(dim(tbl)) - 1 - sum(dim(tbl) - 1), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

pairwise.assoc.test.unconstrained.chisq = function(t0, t1) {
  # The pairwise Pearson chi-squared test on the 2x9 (or less) contingency table.
  # Assuming nothing at all (unconstrained test)
  
  # NOTE: the tables must be pruned so that they do not contain columns with all 
  # zeros.
  
  ret = .my.chisq.test.2d(rbind(as.vector(t0), as.vector(t1)))
  ret$pen = t1 / (t0 + t1)
    
  return (ret)
}

pairwise.assoc.test.unconstrained.gsq = function(t0, t1) {
  # The pairwise G test (GLRT conditioned on the observed marginals of the 2x3
  # contingency table), assuming nothing at all (an unconstrained test)

  ret = .my.g.test.2d(rbind(as.vector(t0), as.vector(t1)))
  ret$pen = t1 / (t0 + t1)
  
  return (ret)
}

pairwise.assoc.test.case.only = function(t1) {
  # Assuming LE, and that disease is rare
  # This is the "case-only" 3x3 test.
  
  # NOTE: Since the test statistic used here can be shown to be independent of the
  # observed marginal distributions, this test ignores any kind of "main effect"
  # (i.e., a marginal association). It tests for a pure interactions above and
  # beyond the main effects. Since the disease is assumed to be rare, the controls
  # carry no signal.
  
  return (.my.g.test.2d(t1))
}

pairwise.assoc.test.ind.3d = function(t0, t1) {
  # Assuming LE
  
  # NOTE:
  #
  # - Here we make use of main effects, and allow controls to play a part. So this 
  # test should be superior to the case-only test when main effects exist (but are
  # maybe too small to be detected by a marginal scan, otherwise we wouldn't need 
  # a pairwise test) and/or if the disease prevalence is not too low so controls
  # are actually useful for detecting a deviation from independence.
  #
  # - This test doesn't assume, or, if you will, use the fact that, even under the
  # alternative, the SNPs should be in LE (the assumption that LE holds regardless
  # of whether or not there is association). Instead this test just checks the
  # case-control independence of y, x1, and x2 (a 2x3x3 contingency table test)
  
  return (.my.g.test.3d(array(c(t0, t1), c(3, 3, 2))))
}

pairwise.assoc.test.pure.unconstrained = function(t0, t1) {
  # The logistic regression based GLRT for full model vs. main effects model

  glm.y = cbind(c(t1), c(t0))
  glm.x = expand.grid(x1 = factor(0:2), x2 = factor(0:2))
  
  anv = anova(glm(glm.y ~ glm.x$x1 + glm.x$x2, family = binomial),
              glm(glm.y ~ glm.x$x1 * glm.x$x2, family = binomial), test = 'Chisq')
  
  statistic = anv$Deviance[2]
  p.value = anv$'Pr(>Chi)'[2]
  
  return (list(pen = t1 / (t0 + t1), statistic = statistic, p.value = p.value))
}

pairwise.assoc.test.kpy = function(t0, t1, prevalence, pen.initial = NULL, pxx.initial = NULL) {
  # Assuming only that disease prevalence is known
  
  txx = t0 + t1
  pxx.null.est = txx / sum(txx)
  
  if (is.null(pen.initial)) {
    pars.initial = c(rep(prevalence, 9), pxx.null.est)
  } else {
    pars.initial = c(pen.initial, pxx.initial)
  }
  
  # NOTE: notice how nice it is that nloptr's interface allows to share code
  # for gradient and objective computations
  
  eval_f = function(pars) {
    pen = pars[1:9]
    pxx = pars[10:18]
    
    objective = -sum(t0 * log(1 - pen) + t1 * log(pen)) - sum(txx * log(pxx))
    gradient = -c(t1 / pen - t0 / (1 - pen), txx / pxx)
    
    return (list(objective = objective, gradient = gradient))
  }
  
  lb = rep(.Machine$double.eps, 18)
  ub = rep(1 - .Machine$double.neg.eps, 18) # not really needed
  
  # TODO: will directly substituting the equality constraint into the target
  # function (but adding an inequality constraint to cover) be faster?
  
  eval_g_eq = function(pars) {
    pen = pars[1:9]
    pxx = pars[10:18]
    
    constraints = c(sum(pen * pxx) - prevalence, sum(pxx) - 1)
    jacobian = matrix(c(pxx, pen, rep(0, 9), rep(1, 9)), nrow = 2, byrow = T)
    
    return (list(constraints = constraints, jacobian = jacobian))
  }
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_GN_ISRES', xtol_abs = 1e-8, maxeval = 100000))
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_LD_AUGLAG', xtol_rel = 1e-8, maxeval = 10000,
  #                          local_opts = list(algorithm = 'NLOPT_LD_LBFGS', xtol_rel = 1e-8)))
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_LD_SLSQP', xtol_rel = 1e-4, xtol_abs = 1e-4 * pars.initial, maxeval = 1000))
  
  sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
               opts = list(algorithm = 'NLOPT_LD_SLSQP', ftol_abs = 1e-4, maxeval = 1000))
  
  pen = matrix(sol$solution[ 1: 9], nrow = 3)
  pxx = matrix(sol$solution[10:18], nrow = 3)
  
  statistic = 2 * (  sum(t0 * log((1 - pen) / (1 - prevalence)) + t1 * log(pen / prevalence), na.rm = T)
                   + sum(txx * log(pxx / pxx.null.est)))
  
  p.value = pchisq(statistic, df = 8, lower.tail = F)
  
  return (list(pen = pen, statistic = statistic, p.value = p.value))
}

pairwise.assoc.test.hwe.le.kpy = function(t0, t1, prevalence, pen.initial = NULL, f1.initial = NULL, f2.initial = NULL) {
  # Assuming HWE and LE, and that disease prevalence is known
  # NOTE: MAFs at both loci have to be estimated.
  
  # Chatterjee and Carrol 2005 (C&C) basically solve the same problem I am in the 
  # first test, but they are interested in the gene-environment problem, where the
  # environmental factor can take many different values. They suggest a
  # semiparametric (or nonparametric) model that optimizes the profile likelihood
  # etc. to make the solution feasible, but this leads to a rather slow
  # implementation in the package "CGEN".
  
  # NOTE: the implementations here assume the contingency tables have SNPs coded
  # 0,1,2 according to count of *minor* allele.
  
  tx1 = rowSums(t0 + t1)
  tx2 = colSums(t0 + t1)
  
  f1.null.est = min(0.5, (0.5 * tx1[2] + tx1[3]) / sum(tx1))
  f2.null.est = min(0.5, (0.5 * tx2[2] + tx2[3]) / sum(tx2))

  if (is.null(pen.initial)) {
    pars.initial = c(rep(prevalence, 9), f1.null.est, f2.null.est)
  } else {
    pars.initial = c(pen.initial, f1.initial, f2.initial)
  }
  
  eval_f = function(pars) {
    pen = pars[1:9]
    f1 = pars[10]
    f2 = pars[11]
    
    p.x1 = c((1 - f1) ^ 2, 2 * f1 * (1 - f1), f1 ^ 2)
    p.x2 = c((1 - f2) ^ 2, 2 * f2 * (1 - f2), f2 ^ 2)

    dp.x1_df1 = c(2 * f1 - 2, 2 - 4 * f1, 2 * f1)
    dp.x2_df2 = c(2 * f2 - 2, 2 - 4 * f2, 2 * f2)
    
    objective = -sum(t0 * log(1 - pen) + t1 * log(pen)) - sum(tx1 * log(p.x1)) - sum(tx2 * log(p.x2))
    gradient = -c(t1 / pen - t0 / (1 - pen), sum(tx1 / p.x1 * dp.x1_df1), sum(tx2 / p.x2 * dp.x2_df2))
    
    return (list(objective = objective, gradient = gradient))
  }

  lb = rep(.Machine$double.eps, 11)
  ub = c(rep(1, 9) - .Machine$double.neg.eps, 0.5, 0.5)
  
  # TODO: will directly substituting the equality constraint into the target
  # function (but adding an inequality constraint to cover) be faster?
  
  eval_g_eq = function(pars) {
    pen = pars[1:9]
    f1 = pars[10]
    f2 = pars[11]
    
    p.x1 = c((1 - f1) ^ 2, 2 * f1 * (1 - f1), f1 ^ 2)
    p.x2 = c((1 - f2) ^ 2, 2 * f2 * (1 - f2), f2 ^ 2)
    
    dp.x1_df1 = c(2 * f1 - 2, 2 - 4 * f1, 2 * f1)
    dp.x2_df2 = c(2 * f2 - 2, 2 - 4 * f2, 2 * f2)
    
    pxx = p.x1 %o% p.x2
    
    constraints = sum(pen * pxx) - prevalence
    jacobian = matrix(c(pxx, sum(pen * (dp.x1_df1 %o% p.x2)), sum(pen * (p.x1 %o% dp.x2_df2))), nrow = 1)
    
    return (list(constraints = constraints, jacobian = jacobian))
  }

# sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
#              opts = list(algorithm = 'NLOPT_LD_AUGLAG', xtol_rel = 1e-4, maxeval = 1000,
#                          local_opts = list(algorithm = 'NLOPT_LD_LBFGS', xtol_rel = 1.0e-4)))
 
  sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
               opts = list(algorithm = 'NLOPT_LD_SLSQP', xtol_abs = 1e-4, maxeval = 1000))
  
  pen = t0 * 0 + sol$solution[1:9]
  f1 = sol$solution[10]
  f2 = sol$solution[11]
  p.x1 = c((1 - f1) ^ 2, 2 * f1 * (1 - f1), f1 ^ 2)
  p.x2 = c((1 - f2) ^ 2, 2 * f2 * (1 - f2), f2 ^ 2)
  
  p1.null.est = c((1 - f1.null.est) ^ 2, 2 * f1.null.est * (1 - f1.null.est), f1.null.est ^ 2)
  p2.null.est = c((1 - f2.null.est) ^ 2, 2 * f2.null.est * (1 - f2.null.est), f2.null.est ^ 2)
  
  statistic = 2 * (  sum(t0 * log((1 - pen) / (1 - prevalence)) + t1 * log(pen / prevalence)) 
                   + sum(tx1 * log(p.x1 / p1.null.est)) + sum(tx2 * log(p.x2 / p2.null.est)))
  
  p.value = pchisq(statistic, df = 8, lower.tail = F)
  
  return (list(pen = pen, statistic = statistic, p.value = p.value))
}

pairwise.assoc.test.kpx.kpy = function(t0, t1, prevalence, pxx, pen.initial = NULL) {
  # Assuming the pairwise distribution, and disease prevalence are known
  
  # We just assume we know the joint distribution of the two SNPs in the 
  # population. This makes sense if we have a lot of prior data (say control 
  # samples from all past studies on the same population), and we assume LE, then
  # genotypic (or allelic + HWE assumption) distributions are easy to estimate.
  
  # NOTE: the notes from the previous test apply here too.
  
  # This is probably still not the most efficient solution
  # I still wonder if I can simply solve this problem directly with some more algebra
  
  if (is.null(pen.initial)) {
    pars.initial = rep(prevalence, 8)
  } else {
    pars.initial = pen.initial[2:9]
  }
  
  fn = function(x) {
    pen = c((prevalence - sum(x * pxx[2:9])) / pxx[1], x)
    return (-sum(t0 * log(1 - pen) + t1 * log(pen)))
  }
  
  gr = function(x) {
    x0 = (prevalence - sum(x * pxx[2:9])) / pxx[1]
    g = t1[2:9] / x - t0[2:9] / (1 - x) - t1[1] / x0 / pxx[1] * pxx[2:9] + t0[1] / (1 - x0) / pxx[1] * pxx[2:9]      
    return (-g)
  }

  ui = rbind(diag(8), -diag(8), -pxx[2:9], pxx[2:9])
  ci = c(rep(.Machine$double.eps, 8), rep(-(1 - .Machine$double.neg.eps), 8), .Machine$double.eps - prevalence, prevalence - pxx[1] * (1 - .Machine$double.neg.eps))
  sol = constrOptim(pars.initial, fn, gr, ui, ci)
  
  pen = t0 * 0 + c((prevalence - sum(sol$par * pxx[2:9])) / pxx[1], sol$par)
  
  statistic = 2 * sum(t0 * log((1 - pen) / (1 - prevalence)) + t1 * log(pen / prevalence))
  p.value = pchisq(statistic, df = length(t0) - 1, lower.tail = F)
  
  return (list(pen = pen, statistic = statistic, p.value = p.value))
}

# "Natural" parameterization. (with nonlinear equality constraint)
pairwise.assoc.test.pop.kpy = function(t0, t1, tp, prevalence, pen.initial = NULL, pxx.initial = NULL) {
  # Assuming only that disease prevalence is known, but we have an extra 
  # population sample.
  
  txx = t0 + t1 + tp
  pxx.null.est = txx / sum(txx)
  null.loglik = sum(t0 * log(1 - prevalence) + t1 * log(prevalence) + txx * log(pxx.null.est), na.rm = T)
  
  if (is.null(pen.initial)) {
    pars.initial = c(rep(prevalence, 9), pxx.null.est)
  } else {
    pars.initial = c(pen.initial, pxx.initial)
  }
  
  # NOTE: notice how nice it is that nloptr's interface allows to share code
  # for gradient and objective computations
  
  eval_f = function(pars) {
    pen = pars[1:9]
    pxx = pars[10:18]
    
    objective = null.loglik - sum(t0 * log(1 - pen) + t1 * log(pen)) - sum(txx * log(pxx))
    gradient = -c(t1 / pen - t0 / (1 - pen), txx / pxx)
    
    return (list(objective = objective, gradient = gradient))
  }
  
  lb = rep(.Machine$double.eps, 18)
  ub = rep(1 - .Machine$double.neg.eps, 18) # not really needed
  
  # TODO: will directly substituting the equality constraint into the target
  # function (but adding an inequality constraint to cover) be faster?
  
  eval_g_eq = function(pars) {
    pen = pars[1:9]
    pxx = pars[10:18]
    
    constraints = c(sum(pen * pxx) - prevalence, sum(pxx) - 1)
    jacobian = matrix(c(pxx, pen, rep(0, 9), rep(1, 9)), nrow = 2, byrow = T)
    
    return (list(constraints = constraints, jacobian = jacobian))
  }
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_GN_ISRES', xtol_abs = 1e-8, maxeval = 100000))
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_LD_AUGLAG', xtol_rel = 1e-8, maxeval = 10000,
  #                          local_opts = list(algorithm = 'NLOPT_LD_LBFGS', xtol_rel = 1e-8)))
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_LD_SLSQP', xtol_rel = 1e-4, xtol_abs = 1e-4 * pars.initial, maxeval = 1000))
  
  #tryCatch(expr={
  sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
               opts = list(algorithm = 'NLOPT_LD_SLSQP', ftol_abs = 1e-8, xtol_rel = 0, maxeval = 1000))
  #}, error = function(e) { browser() })
  
  pen = matrix(sol$solution[ 1: 9], nrow = 3)
  pxx = matrix(sol$solution[10:18], nrow = 3)
  
  statistic = -2 * sol$objective
  p.value = pchisq(statistic, df = 8, lower.tail = F)
  
  return (list(pen = pen, statistic = statistic, p.value = p.value))
}

if (0) {
  # A different parameterization, that has only linear equality constraints
  # I have a problem here, that constrOptim can't deal with the limited domain
  # of the target function. Have I missed something in the box constraints? or
  # is this a real limitation of the barrier method?
  pairwise.assoc.test.pop.kpy = function(t0, t1, tp, prevalence, pen.initial = NULL, pxx.initial = NULL) {
    # Assuming only that disease prevalence is known, but we have an extra 
    # population sample.
    
    txx = t0 + t1 + tp
    pxx.null.est = (txx + 1) / sum(txx + 1)
    pxx.null.est[pxx.null.est == 0] = .Machine$double.eps / prevalence + .Machine$double.eps
    pxx.null.est[pxx.null.est == 1] = 1 - .Machine$double.neg.eps
    qxx.null.est = prevalence * pxx.null.est
    null.loglik = sum(t0 * log(pxx.null.est - qxx.null.est) + t1 * log(qxx.null.est) + tp * log(pxx.null.est))
    
    if (is.null(pen.initial)) {
      pars.initial = c(qxx.null.est[2:9], pxx.null.est[2:9])
    } else {
      pars.initial = c(pen.initial, pxx.initial)
    }
    
    fn = function(pars) {
      q00 = prevalence - sum(pars[1:8])
      p00 = 1 - sum(pars[9:16])
      q = c(q00, pars[1:8])
      p = c(p00, pars[9:16])
      
      null.loglik - sum(t0 * log(p - q) + t1 * log(q) + tp * log(p))
    }
    
    gr = function(pars) {
      q = pars[1:8]
      p = pars[9:16]
      q00 = prevalence - sum(q)
      p00 = 1 - sum(p)
      
      -c( t0[1] / (p00 - q00) - t0[2:9] / (p - q) - t1[1] / q00 + t1[2:9] / q,
         -t0[1] / (p00 - q00) + t0[2:9] / (p - q) - tp[1] / p00 + tp[2:9] / p)
    }
    
    # This is extremely wasteful! but I have no other way to code it
    ui = rbind(diag(16), cbind(-diag(8), diag(8)), c(rep(-1, 8), rep(0, 8)), 
               c(rep(0, 8), rep(-1, 8)), c(rep(1, 8), rep(-1, 8)))
    ci = c(rep(.Machine$double.eps, 24), -prevalence +.Machine$double.eps, 
           -1 + .Machine$double.eps, prevalence - 1 + .Machine$double.eps)
    
    sol = constrOptim(pars.initial, fn, gr, ui, ci)
    
    pxx = matrix(c(1          - sum(sol$par[9:16]), sol$par[9:16]), nrow = 3)
    pen = matrix(c(prevalence - sum(sol$par[1: 8]), sol$par[1: 8]) / pxx, nrow = 3)
    
    statistic = -2 * sol$value
    p.value = pchisq(statistic, df = 8, lower.tail = F)
    
    return (list(pen = pen, statistic = statistic, p.value = p.value))
  }

  # Same parameterization as the previous one, but using NLOpt
  pairwise.assoc.test.pop.kpy = function(t0, t1, tp, prevalence, pen.initial = NULL, pxx.initial = NULL) {
    # Assuming only that disease prevalence is known, but we have an extra 
    # population sample
    
    txx = t0 + t1 + tp
    pxx.null.est = (txx + 1) / sum(txx + 1)
    pxx.null.est[pxx.null.est == 0] = .Machine$double.eps / prevalence + .Machine$double.eps
    pxx.null.est[pxx.null.est == 1] = 1 - .Machine$double.neg.eps
    qxx.null.est = prevalence * pxx.null.est
    null.loglik = sum(t0 * log(pxx.null.est - qxx.null.est) + t1 * log(qxx.null.est) + tp * log(pxx.null.est))
    
    if (is.null(pen.initial)) {
      pars.initial = c(qxx.null.est[2:9], pxx.null.est[2:9])
    } else {
      pars.initial = c(pen.initial, pxx.initial)
    }
    
    eval_f = function(pars) {
      q00 = prevalence - sum(pars[1:8])
      p00 = 1 - sum(pars[9:16])
      q = c(q00, pars[1:8])
      p = c(p00, pars[9:16])
      
      objective = null.loglik - sum(t0 * log(p - q) + t1 * log(q) + tp * log(p))
      gradient = -c( t0[1] / (p[1] - q[1]) - t0[2:9] / (p[2:9] - q[2:9]) - t1[1] / q[1] + t1[2:9] / q[2:9],
                    -t0[1] / (p[1] - q[1]) + t0[2:9] / (p[2:9] - q[2:9]) - tp[1] / p[1] + tp[2:9] / p[2:9])
            
      return (list(objective = objective, gradient = gradient))
    }
    
    lb = rep(.Machine$double.eps, 16)
    ub = rep(Inf, 16) # taken care of by eval_g_ineq (when combined with lb)
        
    eval_g_ineq = function(pars) {
      q00 = prevalence - sum(pars[1:8])
      p00 = 1 - sum(pars[9:16])
      q = c(q00, pars[1:8])
      p = c(p00, pars[9:16])
      sq = sum(q)
      sp = sum(p)
      
      # very wasteful since it is actually a very sparse matrix
      constraints = c(q - p, sq - prevalence, sp - 1, prevalence - 1 - sq + sp)
      jacobian = rbind(c(rep(-1, 8), rep(1, 8)), cbind(diag(8), -diag(8)), 
                       c(rep(1, 8), rep(0, 8)), c(rep(0, 8), rep(1, 8)), c(rep(-1, 8), rep(1, 8)))
      
      return (list(constraints = constraints, jacobian = jacobian))
    }

    sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_ineq = eval_g_ineq, lb = lb, ub = ub,
                 opts = list(algorithm = 'NLOPT_LD_SLSQP', ftol_abs = 1e-8, xtol_rel = 0, maxeval = 1000))
    
    pxx = matrix(c(1          - sum(sol$solution[9:16]), sol$solution[9:16]), nrow = 3)
    pen = matrix(c(prevalence - sum(sol$solution[1: 8]), sol$solution[1: 8]) / pxx, nrow = 3)
    
    statistic = -2 * sol$objective
    p.value = pchisq(statistic, df = 8, lower.tail = F)
    
    return (list(pen = pen, statistic = statistic, p.value = p.value))
  }
}

pairwise.assoc.test.pop.hwe.le.kpy = function(t0, t1, tp1, tp2, prevalence, pen.initial = NULL, f1.initial = NULL, f2.initial = NULL) {
  # Assuming HWE and LE, and that disease prevalence is known. We also
  # have an extra population sample.
  
  tx1 = rowSums(t0 + t1) + tp1
  tx2 = colSums(t0 + t1) + tp2
  
  f1.null.est = min(0.5, (0.5 * tx1[2] + tx1[3]) / sum(tx1))
  f2.null.est = min(0.5, (0.5 * tx2[2] + tx2[3]) / sum(tx2))
  
  p1.null.est = c((1 - f1.null.est) ^ 2, 2 * f1.null.est * (1 - f1.null.est), f1.null.est ^ 2)
  p2.null.est = c((1 - f2.null.est) ^ 2, 2 * f2.null.est * (1 - f2.null.est), f2.null.est ^ 2)
  null.loglik = sum(t0 * log(1 - prevalence) + t1 * log(prevalence)) + sum(tx1 * log(p1.null.est)) + sum(tx2 * log(p2.null.est))
  
  if (is.null(pen.initial)) {
    pars.initial = c(rep(prevalence, 9), f1.null.est, f2.null.est)
  } else {
    pars.initial = c(pen.initial, f1.initial, f2.initial)
  }
  
  eval_f = function(pars) {
    pen = pars[1:9]
    f1 = pars[10]
    f2 = pars[11]
    
    p.x1 = c((1 - f1) ^ 2, 2 * f1 * (1 - f1), f1 ^ 2)
    p.x2 = c((1 - f2) ^ 2, 2 * f2 * (1 - f2), f2 ^ 2)
    
    dp.x1_df1 = c(2 * f1 - 2, 2 - 4 * f1, 2 * f1)
    dp.x2_df2 = c(2 * f2 - 2, 2 - 4 * f2, 2 * f2)
    
    objective = null.loglik - sum(t0 * log(1 - pen) + t1 * log(pen)) - sum(tx1 * log(p.x1)) - sum(tx2 * log(p.x2))
    gradient = -c(t1 / pen - t0 / (1 - pen), sum(tx1 / p.x1 * dp.x1_df1), sum(tx2 / p.x2 * dp.x2_df2))
    
    return (list(objective = objective, gradient = gradient))
  }
  
  lb = rep(.Machine$double.eps, 11)
  ub = c(rep(1, 9) - .Machine$double.neg.eps, 0.5, 0.5)
  
  # TODO: will directly substituting the equality constraint into the target
  # function (but adding an inequality constraint to cover) be faster?
  
  eval_g_eq = function(pars) {
    pen = pars[1:9]
    f1 = pars[10]
    f2 = pars[11]
    
    p.x1 = c((1 - f1) ^ 2, 2 * f1 * (1 - f1), f1 ^ 2)
    p.x2 = c((1 - f2) ^ 2, 2 * f2 * (1 - f2), f2 ^ 2)
    
    dp.x1_df1 = c(2 * f1 - 2, 2 - 4 * f1, 2 * f1)
    dp.x2_df2 = c(2 * f2 - 2, 2 - 4 * f2, 2 * f2)
    
    pxx = p.x1 %o% p.x2
    
    constraints = sum(pen * pxx) - prevalence
    jacobian = matrix(c(pxx, sum(pen * (dp.x1_df1 %o% p.x2)), sum(pen * (p.x1 %o% dp.x2_df2))), nrow = 1)
    
    return (list(constraints = constraints, jacobian = jacobian))
  }
      
# sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
#              opts = list(algorithm = 'NLOPT_GN_ISRES', xtol_abs = 1e-8, maxeval = 100000))
  
# sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
#              opts = list(algorithm = 'NLOPT_LD_AUGLAG', xtol_rel = 1e-8, maxeval = 10000,
#                          local_opts = list(algorithm = 'NLOPT_LD_LBFGS', xtol_rel = 1e-8)))

# sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
#              opts = list(algorithm = 'NLOPT_LD_SLSQP', xtol_rel = 1e-4, xtol_abs = 1e-4 * pars.initial, maxeval = 1000))
  
  sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
               opts = list(algorithm = 'NLOPT_LD_SLSQP', ftol_abs = 1e-8, xtol_rel = 0, maxeval = 1000))
  
  pen = t0 * 0 + sol$solution[1:9]
  f1 = sol$solution[10]
  f2 = sol$solution[11]
  p.x1 = c((1 - f1) ^ 2, 2 * f1 * (1 - f1), f1 ^ 2)
  p.x2 = c((1 - f2) ^ 2, 2 * f2 * (1 - f2), f2 ^ 2)
    
  statistic = -2 * sol$objective
  p.value = pchisq(statistic, df = 8, lower.tail = F)

  return (list(pen = pen, statistic = statistic, p.value = p.value))
}

pairwise.assoc.test.pure.pop.kpy = function(t0, t1, tp, prevalence, pen.initial = NULL, pxx.initial = NULL) {
  # Assuming only that disease prevalence is known, but we have an extra 
  # population sample. This is a test of an unconstrained pairwise association
  # versus a logit-additive model.
  
  txx = t0 + t1 + tp
  
  #
  # Start by fitting the null model
  # logistic MLE y ~ x1 + x2, subject to the known prevalence.
  #
  
  pxx.noassoc.est = txx / sum(txx)
  noassoc.loglik = sum(t0 * log(1 - prevalence) + t1 * log(prevalence) + txx * log(pxx.noassoc.est), na.rm = T)
  
  if (is.null(pxx.initial)) {
    pars.initial = c(.logit(prevalence), rep(0, 4), pxx.noassoc.est)
  } else {
    pars.initial = c(.logit(prevalence), rep(0, 4), pxx.initial) # can I do better?
  }

  A1 = matrix(c(0, 1, 1, 0, 1, 1, 0, 1, 1), 3, 3)
  A2 = matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1), 3, 3)
  B1 = t(A1)
  B2 = t(A2)
  
  eval_f = function(pars) {
    logit_pen = pars[1] + pars[2] * A1 + pars[3] * A2 + pars[4] * B1 + pars[5] * B2
    pen = .logistic(logit_pen)
    pxx = pars[6:14]
    exp_mlogit_pen = exp(-logit_pen)
    
    dl_dpen = t1 / pen - t0 / (1 - pen)
    dl_dlogit_pen = dl_dpen * (1 + exp_mlogit_pen) ^ (-2) * exp_mlogit_pen
    dl_dmu = sum(dl_dlogit_pen)
    dl_da1 = sum(dl_dlogit_pen * A1)
    dl_da2 = sum(dl_dlogit_pen * A2)
    dl_db1 = sum(dl_dlogit_pen * B1)
    dl_db2 = sum(dl_dlogit_pen * B2)
    
    objective = noassoc.loglik - sum(t0 * log(1 - pen) + t1 * log(pen)) - sum(txx * log(pxx))
    gradient = -c(dl_dmu, dl_da1, dl_da2, dl_db1, dl_db2, txx / pxx)
    
    return (list(objective = objective, gradient = gradient))
  }
    
  eval_g_eq = function(pars) {
    logit_pen = pars[1] + pars[2] * A1 + pars[3] * A2 + pars[4] * B1 + pars[5] * B2
    pen = .logistic(logit_pen)
    pxx = pars[6:14]
    exp_mlogit_pen = exp(-logit_pen)
    
    dg_dlogit_pen = pxx * (1 + exp_mlogit_pen) ^ (-2) * exp_mlogit_pen
    dg_dmu = sum(dg_dlogit_pen)
    dg_da1 = sum(dg_dlogit_pen * A1)
    dg_da2 = sum(dg_dlogit_pen * A2)
    dg_db1 = sum(dg_dlogit_pen * B1)
    dg_db2 = sum(dg_dlogit_pen * B2)
    
    constraints = c(sum(pen * pxx) - prevalence, sum(pxx) - 1)
    jacobian = matrix(c(dg_dmu, dg_da1, dg_da2, dg_db1, dg_db2, pen, rep(0, 5), rep(1, 9)), nrow = 2, byrow = T)
    
    return (list(constraints = constraints, jacobian = jacobian))
  }
  
  lb = c(rep(-Inf, 5), rep(.Machine$double.eps, 9))
  ub = rep(Inf, 14)
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_GN_ISRES', xtol_abs = 1e-8, maxeval = 100000))
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_LD_AUGLAG', xtol_rel = 1e-8, maxeval = 10000,
  #                          local_opts = list(algorithm = 'NLOPT_LD_LBFGS', xtol_rel = 1e-8)))
  
  # sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
  #              opts = list(algorithm = 'NLOPT_LD_SLSQP', xtol_rel = 1e-4, xtol_abs = 1e-4 * pars.initial, maxeval = 1000))
  
  #tryCatch(expr={
  sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
               opts = list(algorithm = 'NLOPT_LD_SLSQP', ftol_abs = 1e-8, xtol_rel = 0, maxeval = 1000))
  #}, error = function(e) { browser() })
  
  #logit_pen = sol$solution[1] + sol$solution[2] * A1 + sol$solution[3] * A2 + sol$solution[4] * B1 + sol$solution[5] * B2
  #pen = .logistic(logit_pen)
  #pxx = matrix(sol$solution[6:14], nrow = 3)
  
  altern.sol = pairwise.assoc.test.pop.kpy(t0, t1, tp, prevalence, pen.initial, pxx.initial)
  
  statistic = altern.sol$statistic + 2 * sol$objective
  p.value = pchisq(statistic, df = 4, lower.tail = F)
  
  return (list(pen = altern.sol$pen, statistic = statistic, p.value = p.value))
}

pairwise.assoc.test.pure.pop.hwe.le.kpy = function(t0, t1, tp1, tp2, prevalence, pen.initial = NULL, f1.initial = NULL, f2.initial = NULL) {
  # Assuming HWE and LE, and that disease prevalence is known. We also
  # have an extra population sample. Test for a pure interaction, above and 
  # beyond a (logistic) main-effects model.
  
  tx1 = rowSums(t0 + t1) + tp1
  tx2 = colSums(t0 + t1) + tp2
  
  f1.noassoc.est = min(0.5, (0.5 * tx1[2] + tx1[3]) / sum(tx1))
  f2.noassoc.est = min(0.5, (0.5 * tx2[2] + tx2[3]) / sum(tx2))
  
  p1.noassoc.est = c((1 - f1.noassoc.est) ^ 2, 2 * f1.noassoc.est * (1 - f1.noassoc.est), f1.noassoc.est ^ 2)
  p2.noassoc.est = c((1 - f2.noassoc.est) ^ 2, 2 * f2.noassoc.est * (1 - f2.noassoc.est), f2.noassoc.est ^ 2)
  noassoc.loglik = sum(t0 * log(1 - prevalence) + t1 * log(prevalence)) + sum(tx1 * log(p1.noassoc.est) + tx2 * log(p2.noassoc.est))

  if (is.null(pen.initial)) {
    pars.initial = c(.logit(prevalence), rep(0, 4), f1.noassoc.est, f2.noassoc.est)
  } else {
    pars.initial = c(.logit(prevalence), rep(0, 4), f1.initial, f2.initial) # can I do better?
  }

  A1 = matrix(c(0, 1, 1, 0, 1, 1, 0, 1, 1), 3, 3)
  A2 = matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1), 3, 3)
  B1 = t(A1)
  B2 = t(A2)
  
  eval_f = function(pars) {
    logit_pen = pars[1] + pars[2] * A1 + pars[3] * A2 + pars[4] * B1 + pars[5] * B2
    pen = .logistic(logit_pen)
    f1 = pars[6]
    f2 = pars[7]
    
    p.x1 = c((1 - f1) ^ 2, 2 * f1 * (1 - f1), f1 ^ 2)
    p.x2 = c((1 - f2) ^ 2, 2 * f2 * (1 - f2), f2 ^ 2)
    
    dp.x1_df1 = c(2 * f1 - 2, 2 - 4 * f1, 2 * f1)
    dp.x2_df2 = c(2 * f2 - 2, 2 - 4 * f2, 2 * f2)

    exp_mlogit_pen = exp(-logit_pen)
    dl_dpen = t1 / pen - t0 / (1 - pen)
    dl_dlogit_pen = dl_dpen * (1 + exp_mlogit_pen) ^ (-2) * exp_mlogit_pen
    dl_dmu = sum(dl_dlogit_pen)
    dl_da1 = sum(dl_dlogit_pen * A1)
    dl_da2 = sum(dl_dlogit_pen * A2)
    dl_db1 = sum(dl_dlogit_pen * B1)
    dl_db2 = sum(dl_dlogit_pen * B2)
    
    objective = noassoc.loglik - sum(t0 * log(1 - pen) + t1 * log(pen)) - sum(tx1 * log(p.x1) + tx2 * log(p.x2))
    gradient = -c(dl_dmu, dl_da1, dl_da2, dl_db1, dl_db2, sum(tx1 / p.x1 * dp.x1_df1), sum(tx2 / p.x2 * dp.x2_df2))
    
    return (list(objective = objective, gradient = gradient))
  }
    
  eval_g_eq = function(pars) {
    logit_pen = pars[1] + pars[2] * A1 + pars[3] * A2 + pars[4] * B1 + pars[5] * B2
    pen = .logistic(logit_pen)
    f1 = pars[6]
    f2 = pars[7]
    
    p.x1 = c((1 - f1) ^ 2, 2 * f1 * (1 - f1), f1 ^ 2)
    p.x2 = c((1 - f2) ^ 2, 2 * f2 * (1 - f2), f2 ^ 2)
    pxx = p.x1 %o% p.x2
    
    dp.x1_df1 = c(2 * f1 - 2, 2 - 4 * f1, 2 * f1)
    dp.x2_df2 = c(2 * f2 - 2, 2 - 4 * f2, 2 * f2)

    exp_mlogit_pen = exp(-logit_pen)

    dg_dlogit_pen = pxx * (1 + exp_mlogit_pen) ^ (-2) * exp_mlogit_pen
    dg_dmu = sum(dg_dlogit_pen)
    dg_da1 = sum(dg_dlogit_pen * A1)
    dg_da2 = sum(dg_dlogit_pen * A2)
    dg_db1 = sum(dg_dlogit_pen * B1)
    dg_db2 = sum(dg_dlogit_pen * B2)
        
    constraints = sum(pen * pxx) - prevalence
    jacobian = matrix(c(dg_dmu, dg_da1, dg_da2, dg_db1, dg_db2,
                        sum(pen * (dp.x1_df1 %o% p.x2)), 
                        sum(pen * (p.x1 %o% dp.x2_df2))), nrow = 1)
    
    return (list(constraints = constraints, jacobian = jacobian))
  }
  
  lb = c(rep(-Inf, 5), rep(.Machine$double.eps, 2))
  ub = c(rep(Inf, 5), 0.5, 0.5)
  
  null.sol = nloptr(x0 = pars.initial, eval_f = eval_f, eval_g_eq = eval_g_eq, lb = lb, ub = ub,
                    opts = list(algorithm = 'NLOPT_LD_SLSQP', ftol_abs = 1e-8, xtol_rel = 0, maxeval = 1000))

  altern.sol = pairwise.assoc.test.pop.hwe.le.kpy(t0, t1, tp1, tp2, prevalence, pen.initial, f1.initial, f2.initial)
    
  statistic = altern.sol$statistic + 2 * null.sol$objective
  p.value = pchisq(statistic, df = 4, lower.tail = F)
  
  return (list(pen = altern.sol$pen, statistic = statistic, p.value = p.value))
}

conditional.assoc.test.pure.pop.hwe.le.kpy = function(t0, t1, tp1, tp2, prevalence, pen.initial = NULL, f1.initial = NULL, f2.initial = NULL) {
  # Assuming HWE and LE, and that disease prevalence is known. We also
  # have an extra population sample. Test for a pairwise association, above and 
  # beyond a marginal association of x1 alone, or of x2 alone (the *maximum*
  # p-value of these two options is returned).

  # FIXME: is the null likelihood definition for the marginal and pairwise tests used below the same?
  # FIXME take care of passing the initial conditions to the marginal functions
  
  # 1. Fit the full model using x1 alone
  lambda.x1 = marginal.assoc.test.pop.hwe.kpy(rowSums(t0), rowSums(t1), tp1, prevalence)$statistic

  # 2. Fit the full model using x2 alone
  lambda.x2 = marginal.assoc.test.pop.hwe.kpy(colSums(t0), colSums(t1), tp2, prevalence)$statistic
  
  # 3. Fit the full model using x1 and x2
  lambda.x1x2 = pairwise.assoc.test.pop.hwe.le.kpy(t0, t1, tp1, tp2, prevalence, pen.initial, f1.initial, f2.initial)$statistic
  
  # 4. Perform the GLRT for 1. vs 3. and for 2. vs 3., and return the maximum p-value
  statistic = lambda.x1x2 - max(lambda.x1, lambda.x2)
  p.value = pchisq(statistic, df = 6, lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}
