pscoast <- function(cmd, file=getOption("gmt.file"))
{
  if(is.null(file)) stop("Please pass a valid 'file' argument, or run gmt(file=\"myfile\").")
  owd <- setwd(dirname(file)); on.exit(setwd(owd))

  gmt.system(paste("pscoast",cmd), file=file)

  invisible(NULL)
}
