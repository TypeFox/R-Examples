WriteOut <- function(x, writetype = c("all", "BioGeoBEARS", "coexistence", "graphs", 
                                      "maps", "nexus", "statistics"), areanames = NULL) {
  
  match.arg(writetype)
  
  if (class(x) == "list") {
    if (length(areanames) == 0) {
      areanam <- x[[1]]$areanam
    } else {
      areanam <- areanames
    }
    if (writetype[1] == "all") {
      .NexusOut(x)
      for (i in 1:length(x)) {
        .WriteTablesSpGeo(x[[i]], prefix = names(x)[i])
      }
      if (length(dim(x[[1]]$coexistence_classified)) == 0) {
        warning("no coexistence matrix found")
      } else {
        for (i in 1:length(x)) {
          .OutHeatCoEx(x[[i]], prefix = names(x)[i])
        }
      }
      for (i in 1:length(x)) {
        Spgc2BioGeoBEARS(x[[i]], file = paste(names(x)[i], "_BioGeoBEARS.txt", 
                                         sep = ""))
        .OutPlotSpPoly(x[[i]], prefix = names(x)[i])
        .OutBarChartPoly(x[[i]], prefix = names(x)[i])
        .OutBarChartSpec(x[[i]], prefix = names(x)[i])
        .OutMapAll(x[[i]], prefix = names(x)[i])
        .OutMapPerSpecies(x[[i]], prefix = names(x)[i])
        .OutMapPerPoly(x[[i]], prefix = names(x)[i])
      }
    }
    if (writetype[1] == "graphs") {
      for (i in 1:length(x)) {
        .OutPlotSpPoly(x[[i]], prefix = names(x)[i])
        .OutBarChartPoly(x[[i]], prefix = names(x)[i])
        .OutBarChartSpec(x[[i]], prefix = names(x)[i])
      }
    }
    if (writetype[1] == "maps") {
      for (i in 1:length(x)) {
        .OutMapAll(x[[i]], prefix = names(x)[i])
        .OutMapPerSpecies(x[[i]], prefix = names(x)[i])
        .OutMapPerPoly(x[[i]], areanames = areanames, prefix = names(x)[i])
      }
    }
    if (writetype[1] == "statistics") {
      for (i in 1:length(x)) {
        .WriteTablesSpGeo(x[[i]], prefix = names(x)[i])
      }
    }
    if (writetype[1] == "BioGeoBEARS") {
      for (i in 1:length(x)) {
        Spgc2BioGeoBEARS(x, file = paste(names(x)[i], "BioGeoBEARS.txt", sep = ""))
      }
    }
    if (writetype[1] == "nexus") {
      .NexusOut(x)
    }
    if (writetype[1] == "coexistence") {
      if (length(dim(x[[1]]$coexistence_classified)) == 0) {
        print("No coexistence matrix found")
      } else {
        for (i in 1:length(x)) {
          .OutHeatCoEx(x[[i]], prefix = names(x)[i])
        }
      }
    }
  } else {
    if (length(areanames) == 0) {
      areanam <- x$areanam
    } else {
      areanam <- areanames
    }
    if (writetype[1] == "all") {
      .NexusOut(x)
      .WriteTablesSpGeo(x)
      
      if (length(dim(x$coexistence_classified)) == 0) {
        warning("No coexistence matrix found")
      } else {
        .OutHeatCoEx(x, prefix = "")
      }
      Spgc2BioGeoBEARS(x, file = "BioGeoBEARS.txt")
      .OutPlotSpPoly(x, prefix = "")
      .OutBarChartPoly(x, prefix = "")
      .OutBarChartSpec(x, prefix = "")
      .OutMapAll(x, prefix = "")
      .OutMapPerSpecies(x, prefix = "")
      .OutMapPerPoly(x, areanames = areanam, prefix = "")
    }
    if (writetype[1] == "graphs") {
      .OutPlotSpPoly(x, prefix = "")
      .OutBarChartPoly(x, prefix = "")
      .OutBarChartSpec(x, prefix = "")
      .OutMapAll(x, prefix = "")
      .OutMapPerSpecies(x, prefix = "")
      .OutMapPerPoly(x, areanames = areanam, prefix = "")
    }
    if (writetype[1] == "maps") {
        .OutMapAll(x, prefix = "")
        .OutMapPerSpecies(x, prefix = "")
        .OutMapPerPoly(x, areanames = areanames, prefix = "")
    }
    
    if (writetype[1] == "statistics") {
      .WriteTablesSpGeo(x)
    }
    if (writetype[1] == "nexus") {
      .NexusOut(x)
    }
    if (writetype[1] == "BioGeoBEARS") {
      Spgc2BioGeoBEARS(x, file = "BioGeoBEARS.txt")
    }
    
    if (writetype[1] == "coexistence") {
      if (length(dim(x$coexistence_classified)) == 0) {
        print("No coexistence matrix found")
      } else {
        .OutHeatCoEx(x, prefix = "")
      }
    }
  }
}