modelExpandParam <-
function (model, params) {
#   if (is.GPModel(model))
#     return (modelExpandParam(modelStruct(model), params))

  if ( is.list(params) )
    params <- params$values

  if ( "paramGroups" %in% names(model) )
    params <- params %*% t(model$paramGroups)

  funcName <- paste(model$type, "ExpandParam", sep="")
  func <- get(funcName, mode="function")
  model <- func(model, params)

  return (model)
}
