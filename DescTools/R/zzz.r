.onLoad <- function(libname, pkgname)
{
  options(lastWord = NULL)
  options(lastPP = NULL)

  # packageStartupMessage("Hello all!", " ", domain = "DescTools", appendLF = FALSE)
  invisible()
}
