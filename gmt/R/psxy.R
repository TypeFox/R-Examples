psxy <- function(x, cmd="-J -R -Scp -W2p -O -K", file=getOption("gmt.file"))
{
  if(is.null(file)) stop("Please pass a valid 'file' argument, or run gmt(file=\"myfile\").")
  owd <- setwd(dirname(file)); on.exit(setwd(owd))

  tmp <- paste(dirname(tempdir()), "xy.gmt", sep="/")
  r2gmt(x, tmp)
  gmt.system(paste("psxy",tmp,cmd), file=file, append=TRUE)

  invisible(NULL)
}
