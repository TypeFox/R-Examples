## Convert object name <-> filename, e.g. aBcD <-> dir/a@Bc@D.RData ("@" needed for Windows):
g.data.mash   <- function(dir, obj)
  file.path(dir, paste(gsub("([[:upper:]])", "@\\1", obj), "RData", sep="."))
