#' Import Luminex data
#'
#'
#' @description Import xPONENT software raw data from Luminex assays and 
#' converts in an easy-to-work file
#'
#' @param x a file or directory of files for Fluorescence and Metadata-type or 
#' Bead files respectively
#' @param ... other parameters of the function
#'
#' 
#' @seealso \code{\link{lum_export}} 
#' 
#' @return If import data is Fluorescence-type data several objects are 
#' returned in a list format:
#' \itemize{
#' \item{\code{dtblock}}{, blocks of information from raw data} 
#' \item{\code{raw_metadata}}{, raw metadata}  
#' \item{\code{well_vars}}{, list of variables for the well data}
#' \item{\code{scurve_vars}}{, list of variables for the standard curve data}
#' \item{\code{average_vars}}{, list of variables for the average data}  
#' \item{\code{batch_vars}}{, list of variables for the batch data}
#' \item{\code{name_batch}}{, name of the batch} 
#' \item{\code{type_raw_data}}{, Fluorescence}
#' \item{\code{region_vars}}{, list of variables for the region data}
#' \item{\code{sample_vars}}{, list of variables for the sample data}
#' }
#' If import data is Bead-type data two objects are returned in a list format:
#' \itemize{
#' \item{\code{bead_files}}{, all information from bead csv files}
#' \item{\code{name_batch}}{, name of the batch}
#' \item{\code{type_raw_data}}{, Bead}
#' }
#' 
#' @details The list of variables that return is based on the main structure of 
#' luminex information.
#' 
#' 
#' 
#' @examples
#' imp_path <-  system.file(c("inst","extdata"),"plate1.csv", 
#' package="drLumi")
#' imp <- lum_import(imp_path)
#' imp
#' 
#' @export
lum_import <- function(x, ...) {
    UseMethod("lum_import") 
} 

