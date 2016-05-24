.onLoad <- function(libname, pkgname) {
  # Error if gpg is not installed or cannot be run from the terminal
  if(!system2("gpg", "--version", stdout = FALSE, stderr = FALSE) == 0){
    stop(
      "Please install GPG first or check if GPG can be run from the command line.
      The rcrypt package is just an interface to GPG. See https://gnupg.org/ for
      installation guidelines."
    )
  }
}
