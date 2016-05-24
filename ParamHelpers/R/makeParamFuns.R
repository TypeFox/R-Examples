#' @rdname Param
#' @export
makeNumericParam = function(id, lower = -Inf, upper = Inf, allow.inf = FALSE, default, trafo = NULL,
  requires = NULL, tunable = TRUE) {

  assertString(id)
  assertNumber(lower)
  assertNumber(upper)
  if (!is.null(trafo))
    assertFunction(trafo)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  if (upper < lower)
    stop("No possible value!")
  assertLogical(tunable, len = 1L)
  makeParam(id = id, type = "numeric", len = 1L, lower = lower, upper = upper, allow.inf = allow.inf,
    values = NULL, cnames = NULL, default = default, trafo = trafo,
    requires = requires, tunable = tunable)
}

#' @rdname Param
#' @export
makeNumericVectorParam = function(id, len, lower = -Inf, upper = Inf, allow.inf = FALSE, cnames = NULL,
  default, trafo = NULL, requires = NULL, tunable = TRUE) {

  assertString(id)
  len = asInt(len)
  if (is.numeric(lower) && length(lower) == 1)
    lower = rep(lower, len)
  if (is.numeric(upper) && length(upper) == 1)
    upper = rep(upper, len)
  assertNumeric(lower, min.len = 1L, any.missing = FALSE)
  assertNumeric(upper, min.len = 1L, any.missing = FALSE)
  if (!is.null(cnames))
    assertCharacter(cnames, len = len, any.missing = FALSE)
  if (!is.null(trafo))
    assertFunction(trafo)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  if (any(upper < lower))
    stop("No possible value!")
  assertLogical(tunable, len = 1L)
  makeParam(id = id, type = "numericvector", len = len, lower = lower, upper = upper, allow.inf = allow.inf,
    values = NULL, cnames = cnames, default = default, trafo = trafo,
    requires = requires, tunable = tunable)
}

#' @rdname Param
#' @export
makeIntegerParam = function(id, lower = -Inf, upper = Inf, default, trafo = NULL,
  requires = NULL, tunable = TRUE) {

  assertString(id)
  assertNumber(lower)
  assertNumber(upper)
  if (!is.null(trafo))
    assertFunction(trafo)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  if (upper < lower)
    stop("No possible value!")
  assertLogical(tunable, len = 1L)
  makeParam(id = id, type = "integer", len = 1L, lower = lower, upper = upper,
    values = NULL, cnames = NULL, default = default, trafo = trafo,
    requires = requires, tunable = tunable)
}

#' @rdname Param
#' @export
makeIntegerVectorParam = function(id, len, lower = -Inf, upper = Inf, cnames = NULL,
  default, trafo = NULL, requires = NULL, tunable = TRUE) {

  assertString(id)
  len = asInt(len)
  if (is.numeric(lower) && length(lower) == 1)
    lower = rep(lower, len)
  if (is.numeric(upper) && length(upper) == 1)
    upper = rep(upper, len)
  assertNumeric(lower, min.len = 1L, any.missing = FALSE)
  assertNumeric(upper, min.len = 1L, any.missing = FALSE)
  if (!is.null(cnames))
    assertCharacter(cnames, len = len, any.missing = FALSE)
  if (!is.null(trafo))
    assertFunction(trafo)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  if (any(upper < lower))
    stop("No possible value!")
  assertLogical(tunable, len = 1L)
  makeParam(id = id, type = "integervector", len = len, lower = lower, upper = upper,
    values = NULL, cnames = cnames, default = default, trafo = trafo,
    requires = requires, tunable = tunable)
}

#' @rdname Param
#' @export
makeLogicalParam = function(id, default, requires = NULL, tunable = TRUE) {
  assertString(id)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  values = list(TRUE, FALSE)
  names(values) = c("TRUE", "FALSE")
  assertLogical(tunable, len = 1L)
  makeParam(id = id, type = "logical", len = 1L, lower = NULL, upper = NULL,
    values = values, cnames = NULL, default = default, trafo = NULL,
    requires = requires, tunable = tunable)
}

#' @rdname Param
#' @export
makeLogicalVectorParam = function(id, len, cnames = NULL, default, requires = NULL, tunable = TRUE) {
  assertString(id)
  len = asInt(len)
  if (!is.null(cnames))
    assertCharacter(cnames, len = len, any.missing = FALSE)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  values = list(TRUE, FALSE)
  names(values) = c("TRUE", "FALSE")
  assertLogical(tunable, len = 1L)
  makeParam(id = id, type = "logicalvector", len = len, lower = NULL, upper = NULL,
    values = values, cnames = cnames, default = default, trafo = NULL,
    requires = requires, tunable = tunable)
}

#' @rdname Param
#' @export
makeDiscreteParam = function(id, values, trafo = NULL, default, requires = NULL, tunable = TRUE) {
  assertString(id)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  values = checkValuesForDiscreteParam(id, values)
  assertLogical(tunable, len = 1L)
  makeParam(id = id, type = "discrete", len = 1L, lower = NULL, upper = NULL,
    values = values, cnames = NULL, default = default, trafo = trafo,
    requires = requires, tunable = tunable)
}

#' @rdname Param
#' @export
makeDiscreteVectorParam = function(id, len, values, default, requires = NULL, tunable = TRUE) {
  assertString(id)
  len = asInt(len)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  values = checkValuesForDiscreteParam(id, values)
  assertLogical(tunable, len = 1L)
  makeParam(id = id, type = "discretevector", len = len, lower = NULL, upper = NULL,
    values = values, cnames = NULL, default = default, trafo = NULL,
    requires = requires, tunable = tunable)
}

#' @rdname Param
#' @export
makeFunctionParam = function(id, default = default, requires = NULL) {
  assertString(id)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  makeParam(id = id, type = "function", len = 1L, lower = NULL, upper = NULL,
    values = NULL, cnames = NULL, default = default, trafo = NULL,
    requires = requires, tunable = FALSE)
}

#FIXME: what happens if NA is later used for untyped params? because we might interpret this as
# missing value wrt. dependent params
#' @rdname Param
#' @export
makeUntypedParam = function(id, default, requires = NULL, tunable = TRUE) {
  assertString(id)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  makeParam(id = id, type = "untyped", len = 1L, lower = NULL, upper = NULL,
    values = NULL, cnames = NULL, default = default, trafo = NULL,
    requires = requires, tunable = TRUE)
}

#' @rdname Param
#' @export
makeCharacterParam = function(id, default, requires = NULL) {
  assertString(id)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  makeParam(id = id, type = "character", len = 1L, lower = NULL, upper = NULL,
    values = NULL, cnames = NULL, default = default, trafo = NULL,
    requires = requires, tunable = FALSE)
}

#' @rdname Param
#' @export
makeCharacterVectorParam = function(id, len, cnames = NULL, default, requires = NULL) {
  assertString(id)
  len = asInt(len)
  if (!is.null(requires))
    assert(checkClass(requires, "call"), checkClass(requires, "expression"))
  makeParam(id = id, type = "charactervector", len = len, lower = NULL, upper = NULL,
    values = NULL, cnames = cnames, default = default, trafo = NULL,
    requires = requires, tunable = FALSE)
}

##### small helpers #####

checkValuesForDiscreteParam = function(id, values) {
  if (is.vector(values))
    values = as.list(values)
  assertList(values)

  if (length(values) == 0L)
    stopf("No possible value for discrete parameter %s!", id)

  # check that NA does not occur in values, we use that for "missing state" for dependent params
  # make sure that this works for complex object too, cannot be done with simple is.na
  if (any(sapply(values, isScalarNA)))
    stopf("NA is not allowed as a value for discrete parameter %s.\nParamHelpers uses NA as a special value for dependent parameters.", id)

  n = length(values)
  ns = names(values)
  # if names missing, set all to ""
  if (is.null(ns))
    ns = rep("", n)
  # guess missing names
  for (i in seq_len(n)) {
    v = values[[i]]
    if(is.na(ns[i]) || ns[i] == "") {
      if (is.character(v) || is.numeric(v))
        ns[i] = as.character(v)
    }
  }
  names(values) = ns
  if (!isProperlyNamed(values)) {
    stopf("Not all values for parameter '%s' were named and names could not be guessed!", id)
  }

  # check that NA does not occur in value names, see above
  if ("NA" %in% names(values))
    stopf("NA is not allowed as a value name for discrete parameter %s.\nParamHelpers uses NA as a special value for dependent parameters.", id)

  return(values)
}
