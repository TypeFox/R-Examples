GraphForestPlotForEveryOutlook <-
function(table, 
  binary=TRUE, mean.sd=TRUE,
  higher.is.better=FALSE,
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95,  
  sims=1,
  ...){  
  # Summons a forest plot for each one of 10 outlooks.
  # At each turn, all unpublished studies are assigned the same outlook,
  # and a forest plot is generated complete with summary effect sizes
  # for pubs, unpubs, and pubs with unpubs. 
  # 
  # Args:
  #   table: 

  outlooks <- c("very positive", "positive", "no effect", "negative", "very negative", 
                "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL")
  n <- length(outlooks)

  if(binary == TRUE){
    for (i in 1:n){
      outlook <- outlooks[i]
      CalculateAndGraphForestPlot(table,
        binary=TRUE,
        level=level,
        higher.is.better=higher.is.better,
        vpos=vpos, 
        pos=pos, 
        neg=neg, 
        vneg=vneg,
        sims=sims,
        rustlook=outlook,
        ...)
    }      
  } else {
    for (i in 1:n){
      outlook <- outlooks[i]
      CalculateAndGraphForestPlot(table,
        binary=FALSE,
        mean.sd=mean.sd,
        level=level,
        higher.is.better=higher.is.better,                     
        rustlook=outlook,
        ...)
    }  
  }
}
