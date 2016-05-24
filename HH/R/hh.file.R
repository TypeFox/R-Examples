## options(HHfile.ROOT.DIR="c:/HOME/hh")  ## value recommended in Appendix B

hh.file <- function (file) {
  HH.ROOT <- options()$HHfile.ROOT.DIR
  if (is.null(HH.ROOT)) stop('\nDefine
  options(HHfile.ROOT.DIR="c:/HOME/hh")  ## value recommended in Appendix B
before using the hh.file() function.')
  file.path(HH.ROOT, file)
}

hh.file.DOS <- function (file, displayForCutAndPaste=TRUE) {
  filepath.DOS <- gsub("/", "\\\\", hh.file(file))
  if (displayForCutAndPaste) {
    cat(filepath.DOS, "\n")
    invisible(filepath.DOS)
  }
  else
    filepath.DOS
}
