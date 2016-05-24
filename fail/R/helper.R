asFlag = function(x, default, na.ok = FALSE) {
  if (missing(x)) {
    if (!missing(default))
      return(default)
    stopf("Argument %s is missing", deparse(substitute(x)))
  }
  assertFlag(x, na.ok = na.ok, .var.name = deparse(substitute(x)))
  x
}

asKeys = function(.self, keys, len, default) {
  if (missing(keys)) {
    if (!missing(default))
      return(default)
    stop("Keys are missing")
  }

  if (!is.character(keys)) {
    keys = try(as.character(keys))
    if (is.error(keys))
      stop("Keys must be of type character or be convertible to character")
  }

  if (!missing(len) && length(keys) != len)
      stop("Keys must have length ", len)
  if (anyMissing(keys))
    stop("Keys contain NAs")

  # R variable pattern: "^((\\.[[:alpha:]._]+)|([[:alpha:]]+))[[:alnum:]_.]*$"
  pattern = "^[[:alnum:]._-]+$"
  ok = grepl(pattern, keys)
  if (!all(ok))
    stopf("Key '%s' in illegal format, see help", head(keys[!ok], 1L))
  if (!.self$all.files && any(substr(keys, 1L, 1L) == "."))
    stop("Cannot work with hidden files (files starting with a dot) if 'all.files' is set to TRUE.")

  return(keys)
}

checkPath = function(path) {
  qassert(path, "S1")
  if (!file.exists(path) && !dir.create(path, recursive = TRUE))
    stopf("Could not create directory '%s'", path)
  assertDirectory(path, access = "r")
  path
}

checkExtension = function(extension) {
  qassert(extension, "S1")
  if (grepl("[^[:alnum:]]", extension))
    stop("Extension contains illegal characters: ",
      collapse(strsplit(gsub("[[:alnum:]]", "", extension), ""), " "))
  return(extension)
}

checkCollision = function(keys) {
  dups = duplicated(tolower(keys))
  if (any(dups)) {
    warningf("The following keys result in colliding files on case insensitive file systems: %s",
      collapse(keys[dups]))
  }
  invisible(TRUE)
}

checkCollisionNew = function(new, old) {
  dups = new %nin% old & tolower(new) %in% tolower(old)
  if (any(dups))
    warningf("Keys collide on case insensitive file systems: %s", collapse(new[dups]))
  invisible(TRUE)
}

fn2key = function(.self, fn) {
  return(sub(sprintf("\\.%s$", .self$extension), "", fn))
}

key2fn = function(.self, key) {
  return(file.path(.self$path, sprintf("%s.%s", key, .self$extension)))
}

nkeys = function(.self) {
  length(list.files(.self$path, pattern = sprintf("\\.%s$", .self$extension), ignore.case = TRUE, all.files = .self$all.files))
}
