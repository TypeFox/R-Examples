CompileEffectsToBeAssigned <-
function(table, 
  binary=TRUE, mean.sd=FALSE, 
  higher.is.better=TRUE, 
  vpos=NA, pos=NA, neg=NA, vneg=NA, 
  level=95,
  binary.measure="RR", continuous.measure="SMD", 
  summary.measure="SMD", method="DL"){
  # Assigns both (1) effect sizes based on the CI of the summary effect across published studies, 
  # and (2) effect sizes not based that CI. 
  #
  # Args:
  #   table: the data set
  #   binary: True/False - The data is binary/count. 
  #   higher.is.better: T/R - Higher counts/effect sizes are desired.
  #   mean.sd: T/F - The data is in the form of means and standard deviations. 
  #   vpos: The effect size to assign studies with a "very positive" outlook.
  #   pos:  The effect size to assign studies with a "positive" outlook.
  #   neg:  The effect size to assign studies with a "negative" outlook.
  #   vneg: The effect size to assign studies with a "very negative" outlook.
  #   binary.measure: 
  #   continuous.measure: The 
  #     "SMD" = standardized mean difference (default)
  #     "SMDH" = standardized mean difference w/o assuming equal population variances 
  #              in the two groups 
  #
  # Returns: A list of effect sizes. 

  # Dependencies:
  #   Calls: ExtractPublishedStudies(), SetEffectSizesIndependentOfPub(), SetEffectSizesDependentOnPub()
  
  ## Extract published studies
  pub <- ExtractPublishedStudies(table)
  
  ## Assign effects (not based on CI of summary effect across published studies)
  eff.nopub <- SetEffectSizesIndependentOfPub(binary=binary, 
                higher.is.better=higher.is.better,
                vpos=vpos, pos=pos, neg=neg, vneg=vneg)
  ## Assign effects based on CI of summary effect across published studies)
  eff.onpub <- SetEffectSizesDependentOnPub(table, 
                binary=binary, mean.sd=mean.sd, 
                higher.is.better=higher.is.better, 
                level=level, 
                binary.measure=binary.measure, continuous.measure=continuous.measure,
                summary.measure=summary.measure, method=method)
  effect.sizes.list <- as.list(c(eff.nopub, eff.onpub))
  names(effect.sizes.list) <- c("vpos", "pos", "noef", "neg", "vneg",
    "vposcl", "poscl", "curr", "negcl", "vnegcl")
  return(effect.sizes.list)
}
