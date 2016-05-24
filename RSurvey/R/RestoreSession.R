# This function restores local objects within the current R session.

RestoreSession <- function(path, save.objs, fun.call) {

  if (missing(path)) {
    if (exists("Data")) {
      initial.dir <- getwd()
      if("R" %in% dir(path=initial.dir, full.names=FALSE))
        initial.dir <- file.path(initial.dir, "R")

      path <- tcltk::tkchooseDirectory(initialdir=initial.dir,
                                       title="Choose R Source Directory...")
      path <- tcltk::tclvalue(path)
      if (!nzchar(path))
        return()

      Data("default.dir", path)
    } else {
      path <- file.path(getwd(), "R")
      if (!file.exists(path))
        return()
    }
  }

  if (missing(save.objs))
    save.objs <- NULL

  all.objs <- ls(all.names=FALSE, envir=as.environment(1))
  cur.objs <- all.objs[!all.objs %in% save.objs]

  r.files <- list.files(path, pattern="[.][R]$", full.names=TRUE,
                        recursive=TRUE, ignore.case=TRUE)

  grDevices::graphics.off()
  rm(list=cur.objs, envir=as.environment(1))

  err.msg <- "\n"

  cat("\n")
  for (i in r.files) {
    file.name <- substr(basename(i), 1L, nchar(basename(i)) - 2L)
    if (!file.name %in% save.objs) {
      ans <- try(source(i), silent=TRUE)
      if (inherits(ans, "try-error"))
        err.msg <- paste0(err.msg, i, "\n", ans)
      else
        cat(i, "\n")
    }
  }

  cat(err.msg)

  if (!missing(fun.call))
    eval(parse(text=fun.call))()
}
