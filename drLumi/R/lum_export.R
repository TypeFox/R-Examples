#' Extract and generates data.frame from lum_import object
#'
#' Extract, generates and export different information from a lum_import object
#' 
#' 
#' @usage
#' lum_export(x, y = NULL, path = NULL, backname = "ackgro", 
#'     stanname = "tandar", samplename = "Sample",
#'     dilutionname = "Dilution Factor", batchname = "batch", ...)
#'
#' @param x a \code{lum_import} object with all information of the 
#' xPONENT software imported
#' @param y an optional \code{lum_import} object different from \code{x} 
#' @param path the directory name where all files are exported. The directory 
#' is created according to batch name and system date. Default is \code{NULL}
#' @param backname character vector that defines background
#' @param stanname character vector that defines standard
#' @param samplename character vector that defines the name of the samples 
#' variable
#' @param dilutionname character vector that defines the name of the dilution 
#' factor variable
#' @param batchname character vector that defines the name of the batch field
#'  in CSV
#' @param ... other parameters of the function
#' 
#' @details This function returns a list with all datasets from an object of 
#' class \code{lum_import}. Also is possible to export into several csv files 
#' (for Fluorescence-type data) or zip file (for Bead-type data).
#'
#'
#' @examples
#' # Read all data
#' imp_path <-  system.file(c("inst","extdata"),"plate1.csv", 
#' package="drLumi")
#' imp <- lum_import(imp_path)
#' exp <- lum_export(imp)
#' names(exp)
#' 
#' # See variables
#' imp$well_vars
#' 
#' # Select only 2
#' imp$well_vars <- c("Median", "Net MFI")
#' exp <- lum_export(imp)
#' head(exp$well)
#' 
#' 
#' @import chron 
#' @importFrom gdata rename.vars
#' @importFrom reshape melt.data.frame
#'
#' @export
lum_export <- function(x, y = NULL, path = NULL, backname = "ackgro", 
                stanname = "tandar", samplename = "Sample", 
                dilutionname = "Dilution Factor", 
                batchname = "batch", ...){
    xori <- x
    yori <- y
    if(!inherits(x ,"lum_import")){
        stop(" 'x' object must be a 'lumImport' class")
    } 
    if(!is.null(y)){
        if(!inherits(y ,"lum_import")){
            stop(" 'y' object must be a 'lumImport' class")
        } 
        if(x$type_raw_data==y$type_raw_data){
            stop(" 'x' and 'y' must be 'Bead' and 'Fluorescence' type ")
        } 
    }
    if(!is.null(path)){
        if( !file.exists(path)) stop(" 'path' does not exist")
    }
    both <- FALSE
    if(!is.null(x) & !is.null(y)) both <- TRUE
    ans <- list()
    if(both==TRUE){
        if(x$type_raw_data=="Bead") x <- yori
        if(y$type_raw_data=="Fluorescence") y <- xori
        if(x$name_batch!=y$name_batch){
            warning("The name of the batch in 'x' is different  in 'y' ")
        } 
    } else {
        x <- xori
        y <- xori
    }
    if(x$type_raw_data=="Fluorescence"){
        if(length(x$well_vars)>=1){
            ans$well <- wellsimport(x$dtblock, x$well_vars, x$name_batch)  
        }
        if(length(x$scurve_vars)>=1){
            ans$scurve <- scurveimport(x$dtblock, x$scurve_vars, x$name_batch)
        }
        if(length(x$average_vars)>=1){
            ans$average <- avgimport(x$dtblock, x$average_vars, x$name_batch) 
        }
        if(length(x$batch_vars)>=1){
            ans$batch <- batchimport(x$raw_metadata, x$batch_vars, 
                                    x$name_batch) 
        }
        if(length(x$region_vars)>=1){
            ans$region <- regionimport(x$dtblock, x$region_vars, x$name_batch)
        }
        if(length(x$sample_vars)>=1){
            ans$sample <- sampleimport(x$dtblock, x$sample_vars, x$name_batch, 
                            backname, stanname, samplename, 
                            dilutionname, batchname)    
        }
        ans$name_batch <- x$name_batch
    }
    if(y$type_raw_data=="Bead"){
        ans$bead <- y$bead    
    }
    if(!is.null(path)){
        newpath <- file.path(path, paste(x$name_batch,Sys.Date(), sep="_") )
        dir.create(newpath)  
        name <- names(ans)
        for(i in 1:length(name)){
            write.csv(ans[[i]],
                    file = file.path(newpath,paste(name[i],".csv",sep="")),
                    row.names=FALSE, na = "")
        }
        if(file.exists( file.path(newpath,"bead.csv") ) ){
            zip(file.path(newpath,"bead.zip") , 
            file.path( newpath,"bead.csv"), extras="-j" )
            unlink(file.path( newpath,"bead.csv"))
        }
    }
    class(ans) <- "lum_export"
    ans
}


