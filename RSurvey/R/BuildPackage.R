# Writes a series of commands for creating an R package to a Windows batch file.
# Software requirements for building R packages in a Windows OS include:
#   Rtools;      http://cran.r-project.org/bin/windows/Rtools/
#   MiKTeX;      http://www.miktex.org/download
#   QPDF;        http://qpdf.sourceforge.net/
#   Ghostscript; http://www.ghostscript.com/
# Place QPDF and Ghostscript in the 'Path' environmental variable.

BuildPackage <- function(check.cran=FALSE, no.vignettes=FALSE) {
  if (.Platform$OS.type != "windows")
    stop(call.=FALSE, "This function requires a Windows platform.")

  pkg <- basename(getwd())
  file.name <- pkg

  build.option <- "--resave-data"
  check.option <- ""
  if (check.cran) {
    file.name <- paste0(file.name, "-cran")
    check.option <- paste(check.option, "--as-cran")
  }
  if (no.vignettes) {
    file.name <-  paste0(file.name, "-no_vignettes")
    build.option <- paste(build.option, "--no-build-vignettes")
    check.option <- paste(check.option, "--no-build-vignettes")
  }
  if (!check.cran && !no.vignettes)
    check.option <- paste(check.option, "--no-build-vignettes")

  description <- readLines("DESCRIPTION")
  ver <- strsplit(grep("Version:", description, value=TRUE), " ")[[1]][2]

  path.pkg <- shQuote(getwd())
  path.tmp <- shQuote(file.path("C:", pkg))
  path.git <- shQuote(file.path("C:", pkg, ".git"))
  path.tar <- shQuote(paste0("C:/", pkg, "_", ver, ".tar.gz"))
  path.chk <- shQuote(paste0("C:/", pkg, ".Rcheck"))
  path.cmd <- paste0(R.home(component="bin"), "/Rcmd")
  file.zip <- shQuote(paste0(pkg, "_*"))

  cs <- paste(Sys.getenv("COMSPEC"), "/c")

  cmd <- paste("CD /d", path.pkg)
  cmd <- append(cmd, paste0(cs, " RM -f ", pkg, "*"))
  cmd <- append(cmd, "CD /d C:/")
  cmd <- append(cmd, paste(cs, path.cmd, "REMOVE", pkg))
  cmd <- append(cmd, paste(cs, "RMDIR /S /Q", path.tmp))
  cmd <- append(cmd, paste(cs, "CP -r", path.pkg, shQuote("C:/")))
  cmd <- append(cmd, paste(cs, "RMDIR /S /Q", path.chk))
  cmd <- append(cmd, paste(cs, "RMDIR /S /Q", path.git))
  cmd <- append(cmd, paste(cs, path.cmd, "build", path.tmp, build.option))
  cmd <- append(cmd, paste(cs, path.cmd, "check", check.option, path.tar))
  cmd <- append(cmd, paste(cs, path.cmd, "INSTALL --build", path.tar))
  cmd <- append(cmd, paste(cs, "RMDIR /S /Q", path.tmp))
  cmd <- append(cmd, paste(cs, "MOVE /Y", file.zip, path.pkg))
  cmd <- append(cmd, "pause")

  f <- tcl("tk_getSaveFile", defaultextension=".bat",
           title="Save Batch file As", initialfile=paste0(file.name, ".bat"),
           initialdir=file.path(getwd(), ".."))
  f <- as.character(f)
  if (length(f) == 0)
    return()

  cat(cmd, file=f, sep="\n")
}
