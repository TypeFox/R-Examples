ConvertBinaryToEffectSize <-
function(table, measure="RR"){
  # Converts count data to an effect size (such as log relative risk).  
  #
  # Args:
  #   table: A data set table with binary count data for both the ctrl and expt arms.
  #   measure: which effect size measure to calculate
  #     "RR" = log risk ratio (default)
  #     "OR" = log odds ratio 
  #
  # Returns: Adjusted SMD/Hedges g and its variance under a fixed effects model.
  
  # Dependencies:
  #   Calls: package "metafor" to use function escalc()
  #   Callers: binary.table()

  fixedmodel <- escalc(measure=measure, data=table, 
                  append=TRUE,
                  ai=table$expt.events, n1i=table$expt.n, 
                  ci=table$ctrl.events, n2i=table$ctrl.n)
  return(fixedmodel)
}
