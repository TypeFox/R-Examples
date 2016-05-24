loadRData = function(.self, fn) {
  ee = new.env(parent = emptyenv(), hash = FALSE)
  ns = load(fn, envir = ee)
  if (.self$simplify && length(ns) == 1L)
    return(ee[[ns]])
  return(as.list(ee))
}

saveRData = function(.self, fn, key, value) {
  ee = new.env(parent = emptyenv(), hash = FALSE)
  assign(key, value, envir = ee)
  save(list = key, envir = ee, file = fn)
  return(invisible(key))
}

loadR = function(.self, fn) {
  ee = new.env(parent = .GlobalEnv)
  wrap = if (.self$suppressMessages)
    function(expr) suppressMessages(suppressPackageStartupMessages(expr))
  else
    identity
  wrap(sys.source(fn, ee, chdir = TRUE))
  ns = ls(ee, all.names = TRUE)
  if (.self$simplify && length(ns) == 1L)
    return(ee[[ns]])
  return(as.list(ee))
}

saveR = function(.self, fn, key, value) {
  ee = new.env(parent = emptyenv(), hash = FALSE)
  assign(key, value, envir = ee)
  dump(key, file = fn, envir = ee)
}
