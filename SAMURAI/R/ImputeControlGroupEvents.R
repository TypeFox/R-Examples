ImputeControlGroupEvents <-
function(table, rounded=FALSE){
  # Calculates the summary relative risk across published studies, 
  # then uses that to impute the number of events in the control arms of unpublished studies
  # (assuming the rate of events in control arms of unpubs are the same as the rate across pubs).
  # 
  # Args: 
  #   table:  A data set with event counts and sample sizes for both arms of published studies
  #           and with sample sizes for both arms of unpublished studies. 
  #   rounded: T/F - Imputed counts are rounded to the nearest integer. 
  #
  # Returns:  The same data set with imputed event counts in control arms of unpublished studies.

  pub                 <- ExtractPublishedStudies(table)
  sum.pub.ctrl.n      <- sum(pub$ctrl.n)
  sum.pub.ctrl.events <- sum(pub$ctrl.events)
  pub.ctrl.propn      <- sum.pub.ctrl.events / sum.pub.ctrl.n

  unpub               <- ExtractUnpublishedStudies(table)
  unpub$ctrl.events   <- pub.ctrl.propn * unpub$ctrl.n    

  if(rounded == TRUE){
    unpub$ctrl.events <- round(unpub$ctrl.events)
  }
  out <- rbind(pub, unpub)
  return(out)
}
