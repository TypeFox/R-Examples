import.egvs <-
function(filenames, output.names) {
    if (length(filenames) != length(output.names)) 
        stop("The number of in/output map names do not match!\n\n")
    for (i in 1:length(filenames)) {
        execGRASS("r.in.gdal", input = filenames[i], output = output.names[i], flags = c("e", 
            "o", "overwrite", "quiet"))
        cat(paste(filenames[i], " raster map was sucessfully loaded into GRASS as ", 
            output.names[i], " !\n\n", sep = ""))
    }
}
