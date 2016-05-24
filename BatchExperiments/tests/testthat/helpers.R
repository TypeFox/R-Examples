library(BBmisc)
options(BBmisc.ProgressBar.style = "off")
options(BatchJobs.verbose = FALSE)

conf = getConfig()
conf$default.resources$walltime = 120
conf$default.resources$memory = 256
conf$mail.start = conf$mail.done = conf$mail.error = "none"
conf$raise.warnings = FALSE
conf$max.concurrent.jobs = Inf
setConfig(conf = conf)
rm(conf)

cleanup = function() {
  dir = "unittests-files"
  if (!file.exists(dir))
    return(TRUE)

  i = 1L
  repeat {
    if(unlink(dir, recursive=TRUE) == 0L)
      return(TRUE)
    if (i == 6L)
      return(FALSE)

    i = i + 1L
    Sys.sleep(5)
  }
}

tf = function() {
  file.path("unittests-files", basename(tempfile("unittest")))
}

makeTestRegistry = function(...) {
  fd = tf()
  rd = file.path(fd, "registry")
  dir.create(fd, recursive=TRUE, showWarning=FALSE)
  makeExperimentRegistry(
    id = "unittests",
    seed = 1,
    file.dir = rd,
    work.dir="unittests-files",
    ...
  )
}

in.dir = function(dir, expr) {
  old = setwd(dir)
  on.exit(setwd(old))
  force(expr)
}

stopifnot(cleanup())
