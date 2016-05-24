`ensemble.dummy.variables` <- function(
    xcat=NULL, freq.min=50, most.frequent=5,
    overwrite=TRUE, ...
)
{
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(xcat,"RasterLayer") == F) {stop("parameter xcat is expected to be a RasterLayer")}

# get all categories of the layer
    freqs <- raster::freq(xcat)
    freqs <- freqs[is.na(freqs[,1])==F, ]
    all.categories <- freqs[,1]
    replace.frame <- data.frame(id=all.categories, v=rep(0, length(all.categories)))
    colnames(replace.frame)[2] <- paste(names(xcat), "dummy", sep="")
    freqs <- freqs[order(freqs[,2], decreasing=T),]
    cat(paste("\n", "Frequency of the categories", "\n", sep = ""))
    print(freqs)

# freq.min should at minimum exclude the least frequent category
    least.freq <- min(freqs[,2])
    if (least.freq > freq.min) {freq.min <- least.freq}
   
# get only variables with higher frequencies
   
    freqs <- freqs[freqs[,2] > freq.min, ]
    if (most.frequent < 1) {most.frequent <- length(freqs[, 1])}
    if (most.frequent < length(freqs[, 1])) {
        freqs <- freqs[1:most.frequent,]
    }
    new.categories <- freqs[, 1]
    cat(paste("\n", "dummy layers will be created for the following categories", "\n", sep = ""))
    print(new.categories)

# filename of original layer

#    if (! require(tools)) {stop("tools package not available")}
    filename1 <- raster::filename(xcat)
    extension1 <- paste(".", tools::file_ext(filename1),sep="")

# create the new layers
    for (i in 1:length(new.categories)) {
        extension2 <- paste("_", new.categories[i], ".", tools::file_ext(filename1), sep="")
        filename2 <- gsub(pattern=extension1, replacement=extension2, x=filename1)

# use working file to be able to change names
        extension3 <- paste("_working.", tools::file_ext(filename1), sep="")
        filename3 <- gsub(pattern=extension1, replacement=extension3, x=filename1)
        replace.frame1 <- replace.frame
        replace.frame1[replace.frame1[,1]==new.categories[i], 2] <- 1
        new.name <- paste(names(xcat)[1], "_", new.categories[i], sep="")
        names(replace.frame1)[2] <- new.name
        new.x <- raster::subs(xcat, replace.frame1, by="id", which=new.name, subsWithNA=TRUE, filename=filename3, overwrite=overwrite, ...)
        names(new.x) <- new.name
        raster::writeRaster(new.x, filename=filename2, overwrite=overwrite, ...)
    }

# remove the working file
    cat(paste("\n", "Removing temporary file: ", filename3, "\n", sep = ""))
    file.remove(filename3)
}

