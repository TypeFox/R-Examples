drawmap <- function(data = NULL, map, regionvar = 2, plotvar = 3, ...)
{
  plotmap(map, x = data[, plotvar], id = data[, regionvar], ...)
}

