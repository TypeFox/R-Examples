gmt <- function(par=NULL, file="map.eps", style="s", quiet=TRUE)
{
  owd <- setwd(dirname(file)); on.exit(setwd(owd))

  if(is.character(par))
  {
    gmt.system(paste("gmtdefaults -D",style,sep=""), file=".gmtdefaults4")
    gmt.system(paste("gmtset", par))
  }

  options(gmt.file=file)

  if(!quiet)
  {
    gmtdefaults <- gmt.system("gmtdefaults -L")
    cat(paste(gmtdefaults,collapse="\n"), "\n\n")
    print(options("gmt.file"))
  }

  invisible(options("gmt.file"))
}
