# Utility function.
#
# Generates 'gg-plotable' data.frame.
# @param fn [\code{smoof_function}]\cr
#   Target function.
# @param sequences [\code{list}]\cr
#   List of sequences. One sequence for each parameter.
#   Unified with expand.grid.
# @param par.set [\code{ParamSet}]\cr
#   Parameter set.
# @return [\code{data.frame}]
generateDataframeForGGPlot = function(fn, sequences, par.set) {
  data = do.call(expand.grid, sequences)
  colnames(data) = getParamIds(par.set, with.nr = TRUE, repeated = TRUE)
  data.as.list = dfRowsToList(par.set = par.set, df = data)
  data[["y"]] = sapply(data.as.list, function(data.row) {
    if (violatesConstraints(fn, unlist(data.row))) {
      return(NA)
    }
    return(fn(data.row))
  })
  return(data)
}

# Utility function.
#
# Get actual bound if finite or default value for plotting.
# @param bound [\code{numeric(1)}]\cr
#   Numeric bound.
# @param default [\code{numeric(1)}]\cr
#   Default value. Used if bound is infinite.
# @return [\code{numeric(1)}]
getBounds = function(bound, default) {
  if (any(is.infinite(bound)))
    return(rep(default, length(bound)))
  return(bound)
}

# Check if plotting is possible.
#
# @param x [\code{smoof_function}]\cr
#  Smoof function.
# @return Nothing
checkPlotFunParams = function(x) {
  n.params = getNumberOfParameters(x)
  par.set = getParamSet(x)

  if (n.params > 2L) {
    stopf("Only function with up to 2 parameters can be plotted, but your function has %i", n.params)
  }

  if (isMultiobjective(x)) {
    stopf("Plotting of multiobjective functions not possible.")
  }
}

# Map number of params to the corresponding plot function.
#
# @param x [\code{smoof_function}]\cr
#   Smoof Function.
# @param mapping [\code{list}]\cr
#   Mapping from string to function.
# @return [\code{function}]
getInternalPlotFunction = function(x, mapping) {
  n.params = getNumberOfParameters(x)
  par.set = getParamSet(x)

  if (isNumeric(par.set, include.int = FALSE)) {
    if (n.params == 1L) {
      return(mapping[["1Dnumeric"]])
    } else {
      return(mapping[["2Dnumeric"]])
    }
  } else if (hasDiscrete(par.set) & hasNumeric(par.set, include.int = FALSE)) {
    return(mapping[["2DMixed"]])
  }
  stopf("This type of function cannot be plotted.")
}
