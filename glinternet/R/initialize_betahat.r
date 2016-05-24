initialize_betahat = function(activeSet, activeSetOld, betahatOld, numLevels){

  if (is.null(activeSet)) return(NULL)
  groupSizes = get_group_sizes(activeSet, numLevels)
  if (is.null(activeSetOld)) return(rep(0, sum(groupSizes)+1))

  nVars = sapply(activeSet, function(x) if(is.null(x)) 0 else nrow(x))
  nVarsOld = sapply(activeSetOld, function(x) if(is.null(x)) 0 else nrow(x))
  cumGroupSizes = c(0, cumsum(get_group_sizes(activeSetOld, numLevels))) + 1
  indices = lapply(activeSet, function(x) if (!is.null(x)) c(t(x)) else NULL)
  indicesOld = lapply(activeSetOld, function(x) if (!is.null(x)) c(t(x)) else NULL)

  .Call("R_initialize_beta", double(sum(groupSizes)+1), betahatOld, nVars, nVarsOld, cumGroupSizes, numLevels, indices$cat, indicesOld$cat, indices$cont, indicesOld$cont, indices$catcat, indicesOld$catcat, indices$contcont, indicesOld$contcont, indices$catcont, indicesOld$catcont)
}
