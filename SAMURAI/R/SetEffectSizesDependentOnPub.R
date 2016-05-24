SetEffectSizesDependentOnPub <-
function(table, 
  binary=NA, mean.sd=NA, 
  higher.is.better=NA, 
  level=95, 
  binary.measure="RR", continuous.measure="SMD", 
  summary.measure="SMD", method="DL"){
  # Assigns effects that are based on the CI of the summary effect across published studies.
  #
  # Args:
  #   table: the data set
  #   binary: True/False - The data is binary/count. 
  #   mean.sd: T/F - The data is in the form of means and standard deviations. 
  #   higher.is.better: T/R - Higher counts/effect sizes are desired.
  #   binary.measure:     if(binary == TRUE), 
  #                       which effect size measure do we use for the conversion?
  #     "RR" = relative risk
  #     "OR" = odds ratio 
  #   continuous.measure: if(binary == FALSE && mean.sd == TRUE), 
  #                       which effect size measure do we use for the conversion? 
  #     "SMD" = standardized mean difference (default)
  #     "SMDH" = standardized mean difference w/o assuming equal population variances 
  #              in the two groups 
  #   summary.measure: Which effect size measure do we use for the summary effect?
  #   method: "DL" for the DerSimonian & Laird method (1996) (default)
  #   level: confidence level = 1 - alpha
  # Returns: A list of effect sizes. 

  # Dependencies:
  #   Calls: FillInMissingEffectSizeDefaults(), CalculateSummaryEffect()
  #     ExtractPublishedStudies(), ConvertBinaryToEffectSize(), ConvertMeanSDToSMD(),
  
  ## Extract published studies
  pub <- ExtractPublishedStudies(table)
  
  ## Convert count data to log RR and its variance for each study
  if(binary==TRUE){
    pub <- ConvertBinaryToEffectSize(pub, measure=binary.measure) 
  } else if(mean.sd==TRUE){
    pub <- ConvertMeanSDToSMD(pub, measure=continuous.measure) 
  } 
  
  ## Calculate summary effect across all published studies
  pub.summary <- CalculateSummaryEffect(pub, summary.measure=summary.measure, 
                                        method=method, level=level) 
  
  if(binary == TRUE){
    pub.effect <- pub.summary$exp.m
    pub.effect.lcl <- pub.summary$exp.m.lcl
    pub.effect.ucl <- pub.summary$exp.m.ucl
  } else {
    pub.effect <- pub.summary$m
    pub.effect.lcl <- pub.summary$m.lcl
    pub.effect.ucl <- pub.summary$m.ucl
  }
  
  # Calculate 
  halfdown <- 0.5 * (pub.effect + pub.effect.lcl)
  halfup   <- 0.5 * (pub.effect + pub.effect.ucl)    
  if(higher.is.better == TRUE){
    vposcl <- pub.effect.ucl
    poscl  <- halfup
    negcl  <- halfdown
    vnegcl <- pub.effect.lcl    
  } else {
    vposcl <- pub.effect.lcl
    poscl  <- halfdown
    negcl  <- halfup
    vnegcl <- pub.effect.ucl        
  }
  current <- pub.effect
  effect.sizes.list.pub.ci <- c(vposcl,poscl,current,negcl,vnegcl)
  
  return(effect.sizes.list.pub.ci)
}
