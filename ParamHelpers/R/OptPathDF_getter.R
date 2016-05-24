# return logical index for selection via dob and / or eol
getOptPathDobAndEolIndex = function(op, dob = op$env$dob, eol = op$env$eol) {
  op$env$dob %in% dob & op$env$eol %in% eol
}


#' @export
getOptPathLength.OptPathDF = function(op) {
  nrow(op$env$path)
}

#' @export
getOptPathEl.OptPathDF = function(op, index) {
  index = asInt(index)
  n = getOptPathLength(op)
  if (!(index >= 1 && index <= n))
    stop("Index must be between 1 and ", n, "!")
  e = op$env
  path = e$path
  y = unlist(path[index, op$y.names, drop = FALSE])
  # remove y names from path, only consider x
  path = path[, setdiff(colnames(path), op$y.names), drop = FALSE]
  x = dfRowToList(path, op$par.set, index)
  res = list(x = x, y = y, dob = e$dob[index], eol = e$eol[index])
  # if errmsg there, return it
  if (!is.null(e$error.message))
    res$error.message = e$error.message[index]
  if (!is.null(e$exec.time))
    res$exec.time= e$exec.time[index]
  if (!is.null(e$extra))
    res$extra = e$extra[[index]]
  return(res)
}

#' @export
getOptPathX.OptPathDF = function(op, dob = op$env$dob, eol = op$env$eol) {
  return(as.data.frame(op, include.x = TRUE, include.y = FALSE, include.rest = FALSE, dob = dob, eol = eol))
}

#' @export
getOptPathY.OptPathDF = function(op, names, dob = op$env$dob, eol = op$env$eol, drop = TRUE) {
  if (missing(names))
    names = op$y.names
  else
    c(names, subset = op$y.names)
  assertFlag(drop)
  y = as.matrix(op$env$path[getOptPathDobAndEolIndex(op, dob, eol), names, drop = FALSE])
  if (drop && length(names) == 1L)
    y = as.numeric(y)
  return(y)
}

#' @export
getOptPathDOB.OptPathDF = function(op, dob = op$env$dob, eol = op$env$eol) {
  return(op$env$dob[getOptPathDobAndEolIndex(op, dob, eol)])
}

#' @export
getOptPathEOL.OptPathDF = function(op, dob = op$env$dob, eol = op$env$eol) {
  return(op$env$eol[getOptPathDobAndEolIndex(op, dob, eol)])
}

#' @export
getOptPathErrorMessages.OptPathDF = function(op, dob = op$env$dob, eol = op$env$eol) {
  return(op$env$error.message[getOptPathDobAndEolIndex(op, dob, eol)])
}

#' @export
getOptPathExecTimes.OptPathDF = function(op, dob = op$env$dob, eol = op$env$eol) {
  return(op$env$exec.time[getOptPathDobAndEolIndex(op, dob, eol)])
}

#' @export
getOptPathCol.OptPathDF = function(op, name, dob = op$env$dob, eol = op$env$eol) {
  assertString(name)
  if (getOptPathLength(op) == 0L)
    stopf("Trying to return a col from an empty opt.path")
  if (name %in% colnames(op$env$path))
    return(op$env$path[getOptPathDobAndEolIndex(op, dob, eol), name])
  if (name == "dob")
    return(getOptPathDOB(op, dob, eol))
  if (name == "eol")
    return(getOptPathEOL(op, dob, eol))
  if (name == "exec.time")
    return(getOptPathExecTimes(op, dob, eol))
  if (name == "error.message")
    return(getOptPathErrorMessages(op, dob, eol))
  if (name %in% names(op$env$extra[[1]]))
    return(extractSubList(op$env$extra[getOptPathDobAndEolIndex(op, dob, eol)], name))
  stop("The column you specified is not present in the opt.path.")
}

#' @export
getOptPathCols.OptPathDF = function(op, names, dob = op$env$dob, eol = op$env$eol, row.names = NULL) {
  assertCharacter(names, any.missing = FALSE)
  d = as.data.frame(op, dob = dob, eol = eol, row.names = row.names)
  return(d[, names, drop = FALSE])
}
