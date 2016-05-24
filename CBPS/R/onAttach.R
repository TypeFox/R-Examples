".onAttach" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = pkg))
  title <- packageDescription(pkg, lib.loc = mylib)$Title
  ver <- packageDescription(pkg, lib.loc = mylib)$Version
  author <- packageDescription(pkg, lib.loc = mylib)$Author
  packageStartupMessage(pkg, ": ", title, "\nVersion: ", ver, "\nAuthors: ", author, "\n")
}
