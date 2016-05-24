getAlgorithmFilePath = function(file.dir, id) {
  # fix for case-insensitive file names
  id = gsub("([[:upper:]])", "@\\1", id)
  file.path(file.dir, "algorithms", sprintf("%s.RData", id))
}

getProblemFilePaths = function(file.dir, id) {
  # fix for case-insensitive file names
  id = gsub("([[:upper:]])", "@\\1", id)
  parts = c("static", "dynamic")
  setNames(file.path(file.dir, "problems", sprintf("%s_%s.RData", id, parts)), parts)
}
