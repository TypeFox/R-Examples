.onAttach <- function(...) {
  if (!interactive()) return()
  packageStartupMessage("Welcome to classyfire version 0.1-2.")
}
