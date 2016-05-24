stdz.maps <-
function(map.names, output.names = NULL) {
    if (is.null(output.names)) 
        output.names <- paste("std_", map.names, sep = "")
    if (length(map.names) != length(output.names)) 
        stop("The number of in/output map names for the maps do not match!\n\n")
    for (i in 1:length(map.names)) {
        execGRASS("g.region", raster = map.names[i])
        statistic <- execGRASS("r.univar", map = map.names[i], flags = c("g", "quiet"), intern = TRUE, 
            legacyExec = TRUE)
        map.stdev <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "stddev")], 
            "=")[[1]][[2]])
        map.mean <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "mean")], 
            "=")[[1]][[2]])
        execGRASS("r.mapcalc", expression = paste(output.names[i], "=(", map.names[i], "-", 
            map.mean, ")/", map.stdev, sep = ""), flags = c("overwrite", "quiet"), 
            legacyExec = TRUE)
        cat(paste(map.names[i], " raster map was sucessfully standardized as ", output.names[i], 
            "!\n\n", sep = ""))
    }
}
