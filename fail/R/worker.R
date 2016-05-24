Ls = function(.self, pattern = NULL) {
  keys = fn2key(.self, list.files(.self$path, pattern = sprintf("\\.%s$", .self$extension), ignore.case = TRUE, all.files = .self$all.files))
  if (!is.null(pattern))
    keys = keys[grepl(pattern, keys)]
  return(keys)
}

Get = function(.self, key, use.cache) {
  fn = key2fn(.self, key)
  if (!file.exists(fn))
    stopf("File for key '%s' (%s) not found", key, fn)

  if (use.cache) {
    if (!.self$cache$exists(key))
      .self$cache$put(key, .self$loadFun(.self, fn))
    return(.self$cache$get(key))
  }
  return(.self$loadFun(.self, fn))
}

Put = function(.self, ..., keys, li, use.cache) {
  args = argsAsNamedList(...)
  if (missing(keys))
    keys = names2(args)
  keys = c(asKeys(.self, keys, len = length(args)), asKeys(.self, names2(li)))
  args = c(args, as.list(li))

  if (anyMissing(keys))
    stop("Could not determine all key names from input")
  if (anyDuplicated(keys))
    stop("Duplicated key names")

  checkCollisionNew(keys, Ls(.self))

  if (use.cache)
    mapply(.self$cache$put, key = keys, value = args,
      USE.NAMES = FALSE, SIMPLIFY = FALSE)
  else
    .self$cache$rm(keys)

  mapply(.self$saveFun, fn = key2fn(.self, keys), key = keys, value = args,
    MoreArgs = list(.self = .self), USE.NAMES = FALSE, SIMPLIFY = FALSE)
  invisible(keys)
}

Remove = function(.self, keys) {
  w = function(key) {
    .self$cache$rm(key)
    fn = key2fn(.self, key)
    return(file.exists(fn) && file.remove(fn))
  }
  ok = vlapply(keys, w)
  if (!all(ok))
    warningf("Files not removed: %s", collapse(keys[!ok]))
  return(invisible(ok))
}

Apply = function(.self, FUN, ..., keys, use.cache, simplify, use.names) {
  wrapper = function(.key, ...) {
    res = try(FUN(Get(.self, .key, use.cache = use.cache), ...), silent = TRUE)
    if (is.error(res))
      stopf("Error applying function on key '%s': %s", .key, as.character(res))
    return(res)
  }

  FUN = match.fun(FUN)
  return(sapply(keys, wrapper, ..., USE.NAMES = use.names, simplify = simplify))
}

Mapply = function(.self, FUN, ..., keys, use.cache, moreArgs, simplify, use.names) {
  wrapper = function(.key, ...) {
    res = try(FUN(key = .key, value = Get(.self, .key, use.cache = use.cache), ...), silent = TRUE)
    if (is.error(res))
      stopf("Error applying function on key '%s': %s", .key, as.character(res))
    return(res)
  }
  assertFunction(FUN, args = c("key", "value"))
  return(mapply(wrapper, .key = keys, ..., MoreArgs = moreArgs, USE.NAMES = use.names, SIMPLIFY = simplify))
}

Assign = function(.self, keys, envir, use.cache) {
  w = function(key, envir) {
    x = Get(.self, key, use.cache)
    if (.self$simplify)
      assign(key, x, envir = envir)
    else
      mapply(assign, names(x), x, MoreArgs = list(envir = envir), SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  lapply(keys, w, envir = envir)
  return(invisible(keys))
}

Size = function(.self, keys, unit = "b") {
  size = as.integer(file.info(key2fn(.self, keys))$size)
  setNames(size / UNITCONVERT[unit], keys)
}

Clear = function(.self, keys) {
  return(invisible(.self$cache$rm(keys)))
}

Cached = function(.self) {
  return(.self$cache$keys())
}

AsList = function(.self, keys, use.cache) {
  setNames(lapply(keys, Get, .self = .self, use.cache = use.cache), keys)
}

Info = function(.self) {
  return(.self[c("path", "extension", "use.cache", "simplify")])
}
