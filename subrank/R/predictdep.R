predictdep <- function (knownvalues, dependence, smoothing = c("Uniform", "Beta")) 
{
  smoothing <- match.arg(smoothing)
  NbKnownObs = dim(knownvalues)[1]
  SubSampSize = dim(dependence$cop)[1]
  knownvars = intersect(names(knownvalues), dependence$varnames)
  if (length(knownvars) == length(dependence$varnames)) 
    return(knownvalues)
  knownvalues = knownvalues[knownvars]
  NbKnownDims = length(knownvars)
  rankknownvalues = knownvalues
  UnknwonVars = setdiff(dependence$varnames, knownvars)
  NbUnknwonDims = length(UnknwonVars)
  knowndims = match(knownvars, dependence$varnames)
  UnknwonDims = match(UnknwonVars, dependence$varnames)
  rankpredicted = numeric(NbUnknwonDims * NbKnownObs)
  for (var in knownvars) {
    numvar=pmatch(var,dependence$varnames)
    rankknownvalues[var] = dependence$FdR[[numvar]](unlist(knownvalues[var]))
  }
  rankknownvalues = as.numeric(t(as.matrix(rankknownvalues)))
  epsilon=1/(10*NbKnownObs)
  rankknownvalues=pmax(0,pmin(1-epsilon,rankknownvalues))
  if (smoothing == "Uniform") {
    rankknownvalues = floor(rankknownvalues * SubSampSize)
  }
  else {
    rankknownvalues = rbinom(length(rankknownvalues), SubSampSize - 
                               1, rankknownvalues)
  }
  US = runif(NbKnownObs)
  rankpredicted = .Call("InterTir", as.integer(NbKnownObs), 
                        as.integer(NbKnownDims), as.integer(NbUnknwonDims), as.integer(SubSampSize), 
                        as.double(US), as.double(dependence$cop), as.integer(rankknownvalues), 
                        as.integer(knowndims - 1), as.integer(UnknwonDims - 1)) + 
    1
  if (smoothing == "Uniform") {
    rankpredicted = (rankpredicted + runif(NbKnownObs * NbUnknwonDims) - 
                       1)/SubSampSize
  }
  else {
    rankpredicted = rbeta(NbKnownObs * NbUnknwonDims, rankpredicted, 
                          SubSampSize + 1 - rankpredicted)
  }
  rankpredicted = as.data.frame(matrix(data = rankpredicted, 
                                       ncol = NbUnknwonDims, nrow = NbKnownObs, byrow = TRUE))
  names(rankpredicted) = UnknwonVars
  PredictedValues = rankpredicted
  for (var in UnknwonVars) {
    numvar=pmatch(var,dependence$varnames)
    PredictedValues[var] = dependence$FdRinv[[numvar]](unlist(rankpredicted[var]))
  }
  pred = cbind(knownvalues, PredictedValues)
  return(pred)
}
