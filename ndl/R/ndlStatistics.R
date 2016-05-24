ndlStatistics = function(ndl, ...)
{ 
  pdata=acts2probs(ndl$activationMatrix)
  pmatrix=pdata$p
  predictions=pdata$predicted
  n.outcomes = nlevels(as.factor(ndl$cuesOutcomes$Outcomes))
  frequency=ndl$cuesOutcomes$Frequency
  n.data=sum(frequency)
  df.null = n.outcomes * n.data
  n.predictors = length(ndl$weightMatrix)
  df.model = df.null - n.predictors

  statistics = modelStatistics(
      observed = as.character(ndl$cuesOutcomes$Outcomes),
      predicted = predictions,
      frequency = frequency,
      p.values = pmatrix,
      n.data = n.data,
      n.predictors = n.predictors,
      outcomes = levels(as.factor(ndl$cuesOutcomes$Outcomes)),
      p.normalize = TRUE, 
      cross.tabulation = TRUE,
      p.zero.correction = 1e-10   # this can be specified according to the structure of the data
      )

  statistics <- c(list(n.data = n.data, df.null = df.null, df.model = df.model), statistics)

  class(statistics) <- c("ndlStatistics")
  return(statistics)

}
