stopIfLearnerParams = function(par.set) {
  if(any(sapply(par.set$pars, function(x) inherits(x, "LearnerParameter"))))
    stop("No par.set parameter in 'generateDesign' can be of class 'LearnerParameter'!
      Use basic parameters instead to describe you region of interest!")
}

stopIfFunOrUntypedParams = function(par.set) {
  types = getParamTypes(par.set)
  not.ok = intersect(c("untyped", "function", "character", "charactervector"), types)
  if (length(not.ok) > 0L)
    stopf("Parameters of this type are not allowed, but were found: %s", collapse(not.ok))
}

doBasicGenDesignChecks = function(par.set) {
  assertClass(par.set, "ParamSet")
  if (isEmpty(par.set))
    stop("par.set must not be empty!")
  stopIfLearnerParams(par.set)
  stopIfFunOrUntypedParams(par.set)

  lower = getLower(par.set, with.nr = TRUE)
  upper = getUpper(par.set, with.nr = TRUE)

  if (any(is.infinite(c(lower, upper))))
    stop("Finite box constraints required!")

  return(list(lower = lower, upper = upper))
}


getParamNA = function(par, repeated = FALSE) {
  v = switch(par$type,
    numeric = NA_real_,
    numericvector = NA_real_,
    integer = NA_integer_,
    integervector = NA_integer_,
    discrete = NA_character_,
    discretevector = NA_character_,
    logical = NA,
    logicalvector = NA,
    character = NA_character_,
    charactervector = NA_character_
  )
  if (repeated)
    v = rep(v, par$len)
  return(v)
}

# Create a list with the values from getValues, but respect vector parameters and
# create a single list item for each vector component named paramId.number
getParamSetValues = function(par.set) {
  pids1 = getParamIds(par.set, repeated = TRUE, with.nr = TRUE)
  pids2 = getParamIds(par.set, repeated = TRUE, with.nr = FALSE)
  values1 = getValues(par.set)
  values2 = lapply(pids2, function(x) {names(values1[[x]])})
  setNames(values2, pids1)
}

# Makes sure we keep all levels of a discrete vector and also keep them in order.
fixDesignFactors = function(des, par.set) {
  types.df = getParamTypes(par.set, df.cols = TRUE)
  values = getParamSetValues(par.set)
  for (i in which(types.df == "factor")) {
    des[, i] = factor(des[, i], levels = values[[colnames(des)[i]]])
  }
  des
}
