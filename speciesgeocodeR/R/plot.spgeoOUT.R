plot.spgeoOUT <- function(x, plottype = "summary", plotout = F, mode = c("percent", "total"), moreborders = F, areanames = NULL, 
    ...) {
    if (plottype == "summary") {
      if ((max(x$species_coordinates_in$XCOOR) + 180) - (min(x$species_coordinates_in$XCOOR) + 180) > (max(x$species_coordinates_in$YCOOR) + 
          90) - (min(x$species_coordinates_in$YCOOR) + 90)) {
          layout(matrix(c(1, 2, 1, 2), 2, 2))
          .MapAll(x)
          par(...)
          .PlotSpPoly(x)
          layout(matrix(c(1, 1, 1, 1), 2, 2))
      } else {
          layout(matrix(c(1, 1, 2, 2), 2, 2))
          .MapAll(x)
          par(...)
          .PlotSpPoly(x)
          layout(matrix(1, 1, 1))
      }
    }
    if (plottype == "species") {
        .BarChartSpec(x, mode = mode, plotout = plotout)
    }
    if (plottype == "polygons") {
        .BarChartPoly(x, plotout = plotout)
    }
    if (plottype == "speciesrichness") {
        .PlotSpPoly(x)
    }
    if (plottype == "coexistence") {
        if (length(dim(x$coexistence_classified)) != 0) {
            .HeatPlotCoEx(x, plotout = plotout)
        } else {
            warning("no coexistence matrix found")
        }
    }
    if (plottype == "mapspecies") {
        .MapPerSpecies(x, moreborders = moreborders, plotout = plotout)
    }
    if (plottype == "mappolygons") {
        .MapPerPoly(x, areanames = areanames, plotout = plotout)
    }
    if (plottype == "mapunclassified") {
        .MapUnclassified(x, moreborders = moreborders)
    }
    if (plottype == "mapall") {
        .MapAll(x)
    }
} 
