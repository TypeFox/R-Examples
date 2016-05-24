# http://stackoverflow.com/questions/25014492/geom-bar-pictograms-how-to
source_github <- function(u) {
  # load package
  require(RCurl)

  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script),envir=.GlobalEnv)
}

source_github("https://raw.githubusercontent.com/robertgrant/pictogram/master/pictogram.R")
img <- readPNG(system.file("img", "Rlogo.png", package="png"))
pictogram(icon = img, n = c( 12, 35, 7),
          grouplabels=c("12 R logos","35 R logos","7 R logos"))
