CalculateTauSquared <-
function(table, measure="RR", method="DL", level=95){
  # Given a table, calculate tau^2 for each of the following subsets:
  #   published; unpublished; all (published & unpublished)
  #
  # Args:
  #
  # Returns: 

  ## for testing
  #   table=table1; level=95; measure="RR"; method="DL"
  
  # published
  pub <- table[which(table$outlook == "published"), ] 
  if(nrow(pub) > 0){
    pub.agg <-  CalculateSummaryEffect(table=pub, level=level,
                  summary.measure="SMD", method=method) 
    pub.tau2 <- pub.agg$tau2
  } else { 
    pub.agg <- NA
  }

  # unpublished
  unpub <- table[which(table$outlook != "published"), ]
  if(nrow(unpub) > 0){
    unpub.agg <- CalculateSummaryEffect(unpub, level=level,
                  summary.measure="SMD", method=method)  
    unpub.tau2 <- unpub.agg$tau2
  } else { 
    unpub.agg <- NA 
  }
  
  # all
  all.agg <- CalculateSummaryEffect(table=table, level=level,
              summary.measure="SMD", method=method)
  all.tau2 <- all.agg$tau2
  
  return(list(pub=pub.tau2, unpub=unpub.tau2, all=all.tau2))  
}
