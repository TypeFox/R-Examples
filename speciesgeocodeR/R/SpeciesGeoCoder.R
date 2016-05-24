SpeciesGeoCoder <- function(x, y, coex = FALSE, graphs = TRUE, areanames = "", 
                            occ.thresh = 0, elevation = FALSE, threshold,
                            verbose = FALSE, cleaning = FALSE, ...) {
    if (elevation == T) {
        if (class(x) == "character") {
            coords <- read.table(x, sep = "\t", header = T, fill = T, quote = "")
        }
        if (class(x) == "data.frame") {
            coords <- x
        }
        if (verbose == T) {
            cat("Downloading elevation information.\n")
        }
        coords$ele <- GetElevation(coords)
        
        if (max(threshold) > max(coords$ele, na.rm = T))
        {
          threshold2 <- threshold[which(threshold < max(coords$ele, na.rm = T))]
          warning(sprintf("maximum threshold (%s meter) is higher than maximum elevation in the dataset (%s meter); maximum threshold set to %s meter", 
                          max(threshold), max(coords$ele, na.rm = T), max(threshold2)))
          threshold <- threshold2
        }
      
        threshold <- unique(c(0, threshold, max(coords$ele, na.rm = T) + 1))
        coords$cuts <- cut(coords$ele, breaks = threshold, labels = paste(">", threshold[-length(threshold)], sep = ""))
        tt <- split(coords, coords$cuts)
        
        tt <- lapply(tt, function(x) .adjFormat(x))
        ini <- lapply(tt, function(x) ReadPoints(x, y, cleaning = cleaning, ...))
        outo <- lapply(ini, function(x) SpGeoCodH(x, areanames, occ.thresh = occ.thresh))
        
        names(outo) <- gsub(">", "over_", names(outo))
        names(outo) <- paste(names(outo), "_meters", sep = "")
        for (i in 1:length(outo)) {
            if (coex == T) {
                for (i in 1:length(outo)) {
                  outo[[i]] <- CoExClass(outo[[i]])
                  .OutHeatCoEx(outo[[i]], prefix = names(outo)[i])
                }
            }
            
            .WriteTablesSpGeo(outo[[i]], prefix = names(outo)[i])
        }
        .NexusOut(outo)
        
        if (graphs == T) {
            for (i in 1:length(outo)) {
                .OutPlotSpPoly(outo[[i]], prefix = names(outo)[i]) 
                .OutBarChartPoly(outo[[i]], prefix = names(outo)[i])
                .OutBarChartSpec(outo[[i]], prefix = names(outo)[i])
                .OutMapAll(outo[[i]], prefix = names(outo)[i])
                .OutMapPerSpecies(outo[[i]], prefix = names(outo)[i])
                .OutMapPerPoly(outo[[i]], areanames, prefix = names(outo)[i])
            }
        }
        
    } else {
        ini <- ReadPoints(x, y, cleaning = cleaning, ...)
        outo <- SpGeoCodH(ini, areanames, occ.thresh = occ.thresh)
        
        if (coex == T) {
            outo <- CoExClass(outo)
            .OutHeatCoEx(outo, prefix = "")
        }
        
        .WriteTablesSpGeo(outo)
        .NexusOut(outo)
        
        if (graphs == T) {
            .OutPlotSpPoly(outo, prefix = "")
            .OutBarChartPoly(outo, prefix = "")
            .OutBarChartSpec(outo, prefix = "")
            .OutMapAll(outo, prefix = "")
            .OutMapPerSpecies(outo, prefix = "")
            .OutMapPerPoly(outo, areanames, prefix = "")
        }
    }
} 
