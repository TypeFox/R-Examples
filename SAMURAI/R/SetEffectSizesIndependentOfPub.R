SetEffectSizesIndependentOfPub <-
function(
  binary=NA, higher.is.better=NA, 
  vpos=NA, pos=NA, neg=NA, vneg=NA){
  # Assigns effects that are not based on CI of summary effect across published studies.
  #
  # Args:
  #   binary: True/False - The data is binary/count. 
  #   higher.is.better: T/R - Higher counts/effect sizes are desired.
  #   vpos: The effect size to assign studies with a "very positive" outlook.
  #   pos:  The effect size to assign studies with a "positive" outlook.
  #   neg:  The effect size to assign studies with a "negative" outlook.
  #   vneg: The effect size to assign studies with a "very negative" outlook.
  #
  # Returns: A vector of effect sizes. 

  # Dependencies:
  #   Calls: FillInMissingEffectSizeDefaults()

  # Notes:
  #   This function contains lists of preset default values. 
  
  ## Algorithm defaults
  effects.default.rr <- c(3,2,1,1/2,1/3)
  effects.default.smd <- c(0.8,0.3,0,-0.3,-0.8)
  if(binary==TRUE){
    if(higher.is.better==TRUE){
      effects.default <- effects.default.rr
    } else{
      effects.default <- 1/effects.default.rr
    }  
  } else{
    if(higher.is.better==TRUE){
      effects.default <- effects.default.smd
    } else{
      effects.default <- effects.default.smd*(-1)
    }
  }

  ## User definitions
  effects.user <- c(vpos,pos,NA,neg,vneg)

  # Override algorithm defaults with user definitions
  effect.sizes.list <- FillInMissingEffectSizeDefaults(effects.default,effects.user)  
  return(effect.sizes.list)
}
