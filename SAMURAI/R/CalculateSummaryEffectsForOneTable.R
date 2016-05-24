CalculateSummaryEffectsForOneTable <-
function(table, measure="RR", method="DL", level=95){
  # Given a table, calculate the aggregate estimates for each of the following subsets:
  #   published
  #   unpublished
  #   all (published & unpublished)
  # 
  # Args:
  #   table:  A data set with missing figures already imputed 
  #           and with effect sizes already calculated for each study.
  #
  # Returns a table with 3 rows / summary effects: across pubs, across unpubs, across all.
  #
  # Note: output is formatted for the addpoly() function in the metafor package
  
  # published
  pub <- table[which(table$outlook=="published"),]
  if(nrow(pub)>0){
    pub.agg <- CalculateSummaryEffect(pub, summary.measure="SMD", method=method, level=level) 
  } else { pub.agg <- NA}
  
  # summary effect for unpublished studies only
  unpub <- table[which(table$outlook != "published"),]
  if(nrow(unpub)>0){
    unpub.agg <- CalculateSummaryEffect(unpub, summary.measure="SMD", method=method, level=level) 
  } else { unpub.agg <- NA }

  # overall summary effect
  all.agg <- CalculateSummaryEffect(table, summary.measure="SMD", method=method, level=level) 
  
  out <- as.data.frame(rbind(pub.agg,unpub.agg,all.agg))
  return(out)  
}
