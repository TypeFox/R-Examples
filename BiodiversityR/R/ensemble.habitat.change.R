`ensemble.habitat.change` <- function(
    base.map=file.choose(), 
    other.maps=utils::choose.files(),
    change.folder="ensembles/change",
    RASTER.format="raster", RASTER.datatype="INT1U", RASTER.NAflag=255,
    KML.out=FALSE, KML.folder="kml/change",
    KML.maxpixels=100000, KML.blur=10
)
{
#
    change.function <- function(x, y) {return(10*x + y)}
# output resembles binary output based on current 1/0 + predicted 1/0
# 11 = suitable in current and predicted habitat
# 10 = suitable in current but not suitable in predicted habitat
# 01 = 1 = not suitable in current but suitable in predicted habitat
# 00 = 0 = not suitable in current and predicted habitat

# check whether information likely to be count data
    base.raster <- raster::raster(base.map)
    if(length(base.map[grepl("presence", base.map)]) < 1) {
        cat(paste("Warning: base.map is not in a subfolder 'presence'\n", sep=""))
    }
    if(length(other.maps[grepl("presence", base.map)]) < 1) {
        cat(paste("Warning: base.map is not in a subfolder 'presence'\n", sep=""))
    }
#
    dir.create(change.folder, showWarnings = F)
    if(KML.out==T && KML.folder == "kml/change") {
        dir.create("kml", showWarnings = F)
        dir.create("kml/change", showWarnings = F)
    }
    if(KML.out == T && KML.folder != "kml/change") {dir.create(KML.folder, showWarnings = F)}
    base.raster <- raster::raster(base.map)
    base.name <- names(base.raster)
    raster::setMinMax(base.raster)
    if (raster::maxValue(base.raster) > 1) {
        cat(paste("Warning: base.raster has values larger than 1, hence does not provide presence-absence", sep=""))
    }
    frequencies <- data.frame(array(dim=c(length(other.maps)+1,4)))
    rownames(frequencies)[1] <- paste(base.name, " (base.map)", sep="")
    names(frequencies) <- c("never suitable (00=0)", "always suitable (11)", "no longer suitable (10)", "new habitat (01=1)")
    freq1 <- raster::freq(base.raster)
    freq1 <- freq1[is.na(freq1[,1])==F,]
    if(length(freq1[freq1[,1]==0, 2]) > 0) {frequencies[1, 1] <- freq1[freq1[,1]==0, 2]}
    if(length(freq1[freq1[,1]==1, 2]) > 0) {frequencies[1, 2] <- freq1[freq1[,1]==1, 2]}
    for (i in 1:length(other.maps)) {
        other.raster <- raster::raster(other.maps[i])
        raster.name <- names(other.raster)
        raster::setMinMax(other.raster)
        if (raster::maxValue(base.raster) > 1) {
            cat(paste("Warning: other.maps [", i, "] has values larger than 1, hence does not provide presence-absence", sep=""))
        }
        changes.raster <- other.raster
        filename1 <- paste(change.folder, "/", raster.name, sep="")
        changes.raster <- raster::overlay(x=base.raster, y=other.raster, fun = change.function, na.rm = FALSE)
        names(changes.raster) <- raster.name
        raster::writeRaster(changes.raster, filename=filename1, progress='text', format=RASTER.format, overwrite=TRUE, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#  avoid possible problems with saving of names of the raster layers
        raster::writeRaster(changes.raster, filename="working.grd", overwrite=T)
        working.raster <- raster::raster("working.grd")
        names(working.raster) <- raster.name
        raster::writeRaster(working.raster, filename=filename1, progress='text', format=RASTER.format, overwrite=TRUE, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
        rownames(frequencies)[i+1] <- raster.name
        freq1 <- raster::freq(changes.raster)
        freq1 <- freq1[is.na(freq1[,1])==F,]
        if(length(freq1[freq1[,1]==0, 2]) > 0) {frequencies[i+1, 1] <- freq1[freq1[,1]==0, 2]}
        if(length(freq1[freq1[,1]==11, 2]) > 0) {frequencies[i+1, 2] <- freq1[freq1[,1]==11, 2]}
        if(length(freq1[freq1[,1]==10, 2]) > 0) {frequencies[i+1, 3] <- freq1[freq1[,1]==10, 2]}
        if(length(freq1[freq1[,1]==1, 2]) > 0) {frequencies[i+1, 4] <- freq1[freq1[,1]==1, 2]}
        filename2 <- paste(change.folder, "/", raster.name, sep="")
        if (KML.out == T) {
            filename2 <- paste(KML.folder, "/", raster.name, sep="")
            raster::KML(changes.raster, filename=filename2, col=c("black","blue","red","green"), breaks=c(-1, 0, 1, 10, 11),
                colNA=0, blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE)
        }
    }
    cat(paste("\n", "Codes used in newly created rasters", sep=""))
    cat(paste("\n\t", "Code = 11 indicates that the cell is suitable for the base and the other map", sep=""))
    cat(paste("\n\t", "Code = 00 = 0 indicates that the cell is not suitable for the base and the other map", sep=""))
    cat(paste("\n\t", "Code = 10 indicates lost habitat (suitable for the base map, not suitable for the other map)", sep=""))
    cat(paste("\n\t", "Code = 01 = 1 indicates new habitat (not suitable for the base map, suitable for the other map)", sep=""))
    cat(paste("\n\n", "Frequencies (number of cells) of habitat changes", sep=""))
    cat(paste("\n", "(first row documents habitat for base map, i.e. the 'no change' scenario)", "\n\n", sep=""))
    print(frequencies)
    if (KML.out == T) {
        cat(paste("\n", "Colour coding in KML layer", sep=""))
        cat(paste("\n\t", "Colour = green indicates that the cell is suitable for the base and the other map (Code = 11)", sep=""))
        cat(paste("\n\t", "Colour = black indicates that the cell is not suitable for the base and the other map (Code = 0)", sep=""))
        cat(paste("\n\t", "Colour = red indicates lost habitat (Code = 10; suitable for the base map, not suitable for the other map)", sep=""))
        cat(paste("\n\t", "Colour = blue indicates new habitat (Code = 1; not suitable for the base map, suitable for the other map)", "\n", sep=""))
    }
    return(frequencies)
}
