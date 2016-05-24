.onAttach <- function(lib, pkg)
{

  # This may help with encodings in Mac/Linux
  # Sys.setlocale(locale = "UTF-8")
  # Sys.setlocale(locale = "WINDOWS-1252")

  packageStartupMessage("Sotkanet R tools (C) 2013-2015 rOpenGov\nhttps://github.com/ropengov/sotkanet")

}
