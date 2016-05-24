getExtension =
  #
  # getExtension("foo.R")
  # getExtension("foo.tar.gz")
  #
function(name, multiple = FALSE)
{
  gsub(".*\\.([a-zA-Z0-9]+)$", "\\1", name)
}


otherMIMETypes = c("r" = "text/R-code", # perhaps x-application/r-code
                   "svg" = "image/svg+xml",
                   "json" = "application/json")

globalVariables("mimeTypeExtensions")

guessMIMEType =
  #
  # Determine the MIME type
  #

  # guessMIMEType("foo.txt")
  # guessMIMEType("foo.png")
  # guessMIMEType("foo.jpeg")
  # guessMIMEType("foo.Z")
  # guessMIMEType(c("foo.txt", "foo.png", "foo.jpeg", "foo.Z"))
  #
  # No svg in standard database we constructued, so add in via otherMIMETypes.
  #
function(name, default = NA)
{
  data("mimeTypeExtensions", envir = environment())
  ext = getExtension(name)
  ans = mimeTypeExtensions[tolower(ext)]
  if(any(i <- is.na(ans)) )
    ans[i] = otherMIMETypes[tolower(ext[i])]

  if(any(i <- is.na(ans)))
    ans[i] = default

  structure(ans, names = name)
}
