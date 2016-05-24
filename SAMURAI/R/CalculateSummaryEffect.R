CalculateSummaryEffect <-
function(table, level=95,
  summary.measure="SMD", method="DL"){
  # Compute a random-effects summary effect (cSMD) from a fixed-effects model with Hedges g
  #
  # Args:
  #   table: Table of Hedges g and their variances under a fixed effects model
  #   summary.measure: 
  #     "SMD" = standardized mean difference (default)
  #     "SMDH" = standardized mean difference w/o assuming equal population variances 
  #              in the two groups 
  #   method: "DL" for the DerSimonian & Laird method (1996) (default)
  #   level: confidence level = 1 - alpha
  #
  # Returns: Random effects model meta-analytic summary

  # Dependencies: 
  #   Calls: package 'metafor' to use the function rma()
  
  # Notes: 
  #   The ... in this function acts as a 'garbage collector' for runaway parameters upstream.
  
  # to avoid R CMD CHECK NOTE: "no visible binding for global variable"
  yi <- NULL 
  vi <- NULL 
  
  randmodel <- rma(yi, vi, data=table, measure=summary.measure, method=method, level=level)  
  
  # binary: summary log RR, or summary log OR
  # continuous: effect size SMD
  m <- randmodel$b[1]  
  m.se <- randmodel$se[1]
  m.lcl <- randmodel$ci.lb[1]  
  m.ucl <- randmodel$ci.ub[1]
  
  # binary: summary RR, or summary OR
  # continuous: effect size SMD
  expm <- exp(m) 
  expm.lcl <- exp(randmodel$ci.lb[1]) 
  expm.ucl <- exp(randmodel$ci.ub[1]) 
  
  # measures of heterogeneity
  tau2  <- randmodel$tau2[1]
  Q     <- randmodel$QE[1]  
  Qpval <- randmodel$QEp[1]

  out <- as.list(c(m, m.se, m.lcl, m.ucl, expm.lcl, expm, expm.ucl,tau2, Q, Qpval))
  names(out) <- c("m","m.se","m.lcl","m.ucl", "exp.m.lcl","exp.m","exp.m.ucl","tau2","Q","Qpval")    
  return(out)
  # return(list(m=m, m.se=m.se, m.lcl=m.lcl, m.ucl=m.ucl, 
  #             expm.lcl=expm.lcl, expm=expm, expm.ucl=expm.ucl,
  #             tau2=tau2, Q=Q, Qpval=Qpval))
}
