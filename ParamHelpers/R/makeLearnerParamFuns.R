#' @rdname LearnerParam
#' @export
makeNumericLearnerParam = function(id, lower = -Inf, upper = Inf, allow.inf = FALSE, default,
  when = "train", requires = NULL, tunable = TRUE) {

  p = makeNumericParam(id, lower, upper, allow.inf = allow.inf, default = default, requires = requires, tunable = tunable)
  learnerParamFromParam(p, when)
}

#' @rdname LearnerParam
#' @export
makeNumericVectorLearnerParam = function(id, len = as.integer(NA), lower = -Inf,
  upper = Inf, allow.inf = FALSE, default, when = "train", requires = NULL, tunable = TRUE) {

  len = asInt(len, na.ok = TRUE)
  if (is.na(len))
    p = makeVectorParamNALength(makeNumericVectorParam, default = default,
      id = id, lower = lower, upper = upper, allow.inf = allow.inf, requires = requires, tunable = tunable)
  else
    p = makeNumericVectorParam(id, len = len, lower = lower, upper = upper, allow.inf = allow.inf, default =  default,
      requires = requires, tunable = tunable)
  p = learnerParamFromParam(p, when)
  p$len = len
  return(p)
}


#' @rdname LearnerParam
#' @export
makeIntegerLearnerParam = function(id, lower = -Inf, upper = Inf,
  default, when = "train", requires = NULL, tunable = TRUE) {

  p = makeIntegerParam(id, lower, upper, default = default, requires = requires, tunable = tunable)
  learnerParamFromParam(p, when)
}

#' @rdname LearnerParam
#' @export
makeIntegerVectorLearnerParam = function(id, len = as.integer(NA), lower = -Inf,
  upper = Inf, default, when = "train", requires = NULL, tunable = TRUE) {

  len = asInt(len, na.ok = TRUE)
  if (is.na(len))
    p = makeVectorParamNALength(makeIntegerVectorParam, default = default,
      id = id, lower = lower, upper = upper, requires = requires, tunable = tunable)
  else
    p = makeIntegerVectorParam(id, len = len, lower = lower, upper = upper, default = default,
      requires = requires, tunable = tunable)
  p = learnerParamFromParam(p, when)
  p$len = len
  return(p)
}

#' @rdname LearnerParam
#' @export
makeDiscreteLearnerParam = function(id, values, default,
  when = "train", requires = NULL, tunable = TRUE) {

  p = makeDiscreteParam(id, values, default = default, requires = requires, tunable = tunable)
  learnerParamFromParam(p, when)
}

#' @rdname LearnerParam
#' @export
makeDiscreteVectorLearnerParam = function(id, len = as.integer(NA), values, default,
  when = "train", requires = NULL, tunable = TRUE) {

  len = asInt(len, na.ok = TRUE)
  if (is.na(len))
    p = makeVectorParamNALength(makeDiscreteVectorParam, default = default,
      id = id, values = values, requires = requires, tunable = tunable)
  else
    p = makeDiscreteVectorParam(id, len = len, values = values, default = default, requires = requires,
      tunable = tunable)
  p = learnerParamFromParam(p, when)
  p$len = len
  return(p)
}

#' @rdname LearnerParam
#' @export
makeLogicalLearnerParam = function(id, default, when = "train", requires = NULL, tunable = TRUE) {

  p = makeLogicalParam(id, default = default, requires = requires, tunable = tunable)
  learnerParamFromParam(p, when)
}

#' @rdname LearnerParam
#' @export
makeLogicalVectorLearnerParam = function(id, len = as.integer(NA), default, when = "train",
  requires = NULL, tunable = TRUE) {

  len = asInt(len, na.ok = TRUE)
  if (is.na(len))
    p = makeVectorParamNALength(makeLogicalVectorParam, default = default,
      id = id, requires = requires, tunable = tunable)
  else
    p = makeLogicalVectorParam(id, len = len, default = default, requires = requires, tunable = tunable)
  p = learnerParamFromParam(p, when)
  p$len = len
  return(p)
}

#' @rdname LearnerParam
#' @export
makeUntypedLearnerParam = function(id, default, when = "train", requires = NULL, tunable = TRUE) {
  p = makeUntypedParam(id, default = default, requires = requires, tunable = tunable)
  learnerParamFromParam(p, when)
}

#' @rdname LearnerParam
#' @export
makeFunctionLearnerParam = function(id, default, when = "train", requires = NULL) {
  p = makeFunctionParam(id, default = default, requires = requires)
  learnerParamFromParam(p, when)
}

learnerParamFromParam = function(p, when) {
  assertChoice(when, c("train", "predict", "both"))
  makeLearnerParam(p, when)
}

makeVectorParamNALength = function(fun, default, ...)  {
  len = if (missing(default)) 1L else length(default)
  fun(len = len, default = default, ...)
}

