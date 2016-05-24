stitch_dummies <- function(scriptPath, htmlPath){
  require(knitr)
  if(!file.exists(scriptPath)) stop("script path doesn't exist")
  if(!file.exists(htmlPath)) stop("html path doesn't exist")
  scripts <- normalizePath(list.files(scriptPath, pattern=".r$", full.names = TRUE))
  wd <- getwd()
  on.exit(setwd(wd))
  setwd(htmlPath)
  for (scr in scripts){
    message(scr)
    stitch_rhtml(scr)
  }
}

stitch_dummies(scriptPath = "rfordummies/inst/scripts/2-clean", htmlPath="rfordummies/inst/scripts/3-html")
