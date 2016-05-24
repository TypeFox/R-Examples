ConvertMeanSDToSMD <-
function(table, measure="SMD"){
  # Converts mean/standard deviation (sd) data to standardized mean differences (SMD).  
  #
  # Args:
  #   table: A data set table with continuous data (mean,sd,n) for both ctrl and expt arms.
  #   measure: which effect size measure to calculate
  #     "SMD" = standardized mean difference (default)
  #     "SMDH" = standardized mean difference w/o assuming equal population variances 
  #              in the two groups 
  #
  # Returns: Adjusted SMD/Hedges g (yi) and its variance (vi) under a fixed effects model.

  # Dependencies: 
  #   Calls: package "metafor" to use function escalc()
  #   Callers: contintuousforest(), binary.table()

  fixedmodel <- escalc(measure=measure, data=table, 
                  append=TRUE,
                  m1i=table$expt.mean, sd1i=table$expt.sd, n1i=table$expt.n,
                  m2i=table$ctrl.mean, sd2i=table$ctrl.sd, n2i=table$ctrl.n)  
  return(fixedmodel)
}
