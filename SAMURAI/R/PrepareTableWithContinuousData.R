PrepareTableWithContinuousData <-
function(table,
  mean.sd=TRUE, 
  higher.is.better=TRUE, 
  rustlook=NA, 
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95, 
  continuous.measure="SMD", summary.measure="SMD", method="DL", 
  seed=NA, noise=0.01){ 
  # Input a table and impute Hedges' g and its variance for each study.
  # 
  # Args:
  #
  # Returns: A table with all events imputed, and with Hedges' g and its variance calculated.

  ## testing
  #   table <- greentea
  #   binary=F; continuous.measure="SMD"; mean.sd=TRUE; higher.is.better=TRUE;
  #   vpos=NA; pos=NA; neg=NA; vneg=NA; level=95; method="DL"
  #   rustlook="negative"
  
  if(mean.sd==TRUE){
    table <- ConvertMeanSDToSMD(table)  
  }
  
  ## Assign effects (not based on CI of summary effect across published studies)
  effect.sizes.list <- CompileEffectsToBeAssigned(table=table,
                          binary=FALSE,  # since we have continuous data
                          higher.is.better=higher.is.better,
                          vpos=vpos, pos=pos,neg=neg,vneg=vneg,
                          mean.sd=mean.sd, 
                          binary.measure=NA, continuous.measure=continuous.measure,
                          summary.measure=summary.measure,
                          method=method, level=level)                                        
  
  ## If desired, assign all unpublished studies the same outlook
  if(!is.na(rustlook)){
    table <- AssignSameRustlook(table,rustlook)
#     levels(table$outlook)
  }  

  table <- ImputeSMD(table, effect.sizes.list,
                     seed=seed, noise=noise)
  table <- ImputeSMDVariance(table)
  
  return(table)  
}
