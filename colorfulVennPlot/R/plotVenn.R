plotVenn <-
function(x, ...)
{ # chooses which plotVenn function to call based on length of x

  if (missing(x)) plotVenn4d(...)
  else
  {
    if (length(x) == 3) plotVenn2d(x, ...)
    if (length(x) == 7) plotVenn3d(x, ...)
    if (length(x) == 15) plotVenn4d(x, ...)
    if (!length(x) %in% c(3,7,15))
      stop('Specified data cannot be plotted by any of plotVenn diagrams. Please, check your data. Make sure to remove the "000" value.')
  }

}
