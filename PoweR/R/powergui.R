power.gui <- function () {
  path <- path.package(package = "PoweR")
  path.gui <- paste(path, "power.gui.R", sep = .Platform$file.sep)
  source(path.gui)
}
