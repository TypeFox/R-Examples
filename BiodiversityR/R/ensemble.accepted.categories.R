`ensemble.accepted.categories` <- function(
    xcat=NULL, categories=NULL,
    filename=NULL, overwrite=TRUE, ...
)
{
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(xcat,"RasterLayer") == F) {stop("parameter xcat is expected to be a RasterLayer")}
    if(is.null(categories) == T) {stop("accepted categories are missing")}
    if(is.null(filename) == T) {
        cat(paste("\n", "No new filename was provided", sep = ""))
#        if (! require(tools)) {stop("tools package not available")}
        filename1 <- filename(xcat)
        extension1 <- paste(".", tools::file_ext(filename1), sep="")
        extension2 <- paste("_new.", tools::file_ext(filename1), sep="")
        filename <- gsub(pattern=extension1, replacement=extension2, x=filename1)
        cat(paste("\n", "New raster will be saved as: ", filename, "\n", sep = ""))
    }
#
    all.categories <- raster::freq(xcat)[,1]
    all.categories <- all.categories[is.na(all.categories) == F]
    new.categories <- all.categories[is.na(match(all.categories, categories) )]
    cat(paste("\n", "categories that will be reclassified as 'NA'", "\n", sep = ""))
    print(new.categories)
#
    replace.frame <- data.frame(id=categories, v=categories)
    colnames(replace.frame)[2] <- names(xcat)
    new.x <- raster::subs(xcat, replace.frame, by=1, which=2, subsWithNA=TRUE, filename=filename, overwrite=overwrite, ...)
    return(filename)
}
