CalculateAndGraphForestPlot <-
function(table,
  binary=TRUE, mean.sd=FALSE,
  higher.is.better=TRUE, 
  rustlook=NA, vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95,
  binary.measure="RR", continuous.measure="SMD", 
  summary.measure="SMD", method="DL",
  seed=NA, noise=0.01, sims=1, 
  title=NA, digits=3, ...){
  # Given a table with binary event data, 
  # this function imputes data, calculates summary effect sizes, 
  # then graphs a forest plot.
  #
  # Args:
  #   level: confidence level = 1-alpha
  #   vpos: "very positive" outlook
  #
  # Returns: a forest plot for a table of pub & unpub studies with binary outcomes.
  #
  # Notes:
  #   Unlike the function GraphBinaryForestPlot(), this function allows
  #   the user to tweak the default effect sizes to be assigned. 
  
  # Dependencies:
  #   Callers: forestsens()
  
  if(binary == TRUE){
    table1 <- PrepareTableWithBinaryData(table,
                higher.is.better=higher.is.better,
                rustlook=rustlook,
                vpos=vpos, pos=pos, neg=neg, vneg=vneg,
                level=level, binary.measure=binary.measure, 
                summary.measure=summary.measure, method=method, 
                seed=seed, sims=sims)

    GraphBinaryForestPlot(table=table1, 
      level=level, higher.is.better=higher.is.better, 
      digits=digits, title=title, ...)     

  } else {
    table1 <- PrepareTableWithContinuousData(table,
                mean.sd=mean.sd, 
                higher.is.better=higher.is.better, 
                vpos=vpos, pos=pos, neg=neg, vneg=vneg,
                rustlook=rustlook, 
                level=level, continuous.measure=continuous.measure, 
                summary.measure=summary.measure, method=method, 
                seed=seed, noise=noise)

    GraphContinuousForestPlot(table=table1, 
      level=level, higher.is.better=higher.is.better, 
      digits=digits, title=title, ...)     
  }
}
