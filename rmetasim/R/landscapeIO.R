
landscape.clean <- function(rland)
  {
    if (is.landscape(rland))
      .Call("clean_landscape",rland,PACKAGE = "rmetasim")
  }


