.onLoad <- function(libname,pkgname)
{
  options(Hyde_fitModel=FALSE)
  options(Hyde_maxDigits = 5)
  options(Hyde_plotOptions = 
    data.frame(type = c("variable", "determ", "decision", "utility"),
               fillcolor = c("white", "white", "#6BAED6", "#FFFFB2"),
               shape = c("ellipse", "ellipse", "rect", "diamond"),
               fontcolor = c("black", "gray70", "black", "black"),
               color = c("black", "gray70", "black", "black"),
               style = c("filled", "filled", "filled", "filled"),
               stringsAsFactors=FALSE))
}

.onUnload <- function(libPath)
{
  options(Hyde_fitModel=NULL)
  options(Hyde_maxDigits = NULL)
  options(Hyde_plotOptions = NULL)
}