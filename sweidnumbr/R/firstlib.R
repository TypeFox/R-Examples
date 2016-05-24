.onAttach <- function(lib, pkg)
{

  # This may help with encodings in Mac/Linux
  # Sys.setlocale(locale = "UTF-8")
  # Sys.setlocale(locale = "WINDOWS-1252")

  packageStartupMessage("sweidnumbr: R tools to handle swedish identity numbers.\nhttps://github.com/rOpenGov/sweidnumbr\n")

}
