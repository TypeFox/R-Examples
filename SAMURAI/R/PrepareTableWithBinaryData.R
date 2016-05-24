PrepareTableWithBinaryData <-
function(table,
  higher.is.better=NA, 
  rustlook=NA,
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95, 
  binary.measure="RR", summary.measure="SMD", method="DL", 
  seed=NA, sims=1){
  # Input a table and impute all events for each study, then calculate Hedges' g and its variance.
  # 
  # Args:
  #
  # Returns: A table with all events imputed, and with Hedges' g and its variance calculated.

  ## testing
  #   binary=TRUE; binary.measure="RR"; continuous.measure="SMD"; mean.sd=FALSE; higher.is.better=TRUE;
  #   vpos=NA; pos=NA; neg=NA; vneg=NA; rustlook=NA; level=95; method="DL"
  
  ## Avoid R CMD CHECK NOTE: "no visible binding for global variable".
  expt.events <- expt.n <- ctrl.events <- ctrl.n <- NULL    
  
  ## Assign effects (not based on CI of summary effect across published studies)
  effect.sizes.list <- CompileEffectsToBeAssigned(table=table,
                        binary=TRUE, 
                        higher.is.better=higher.is.better,
                        vpos=vpos, pos=pos, neg=neg, vneg=vneg,
                        binary.measure=binary.measure, 
                        summary.measure=summary.measure, method=method, 
                        level=level)
    
  ## If desired, assign all unpublished studies the same outlook
  if(!is.na(rustlook)){
    table <- AssignSameRustlook(table,rustlook)
    #     levels(table$outlook)
  }  

  table <- ImputeBinaryEvents(table,effect.sizes.list,sims=sims)
  table <- ConvertBinaryToEffectSize(table,measure=binary.measure)

  return(table)  
}
