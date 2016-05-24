info = function(...) {
  if (getOption("BatchJobs.verbose", default = TRUE))
    message(sprintf(...))
}

# check for valid algorithm and problem ids
# we can be more relaxed her, this does not affect
# sql table names (unlike the registry id)
checkIdValid = function(id, allow.minus = TRUE) {
  assertString(id)
  if (allow.minus)
    pattern = "^[a-zA-Z]+[0-9a-zA-Z_.-]*$"
  else
    pattern = "^[a-zA-Z]+[0-9a-zA-Z_.]*$"
  if (!grepl(pattern, id))
    stopf("Id does not comply with pattern %s: %s", pattern, id)
}
