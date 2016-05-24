#' @name validate
#' @aliases validate_filename_output
#' @export
#' 
#' @title Functions that check the validty of values throughout the workflow.
#'  
#' @description These functions help identify mistakes in formatting before the create difficult-to-diagnose problems later.
#' 
#' @param filename The name of the file to be validated.
#' @param path The location of the file to be validated.
#' @param file_extension_expected The extension of the file.  This defaults to "out", which corresponds to Mplus output.
#' @param underscore_count_expected The number of underscores required in the name (not currently used).
#' 
#' @return An \code{invisible} \code{TRUE} value if the filename is valid.  Otherwise, an error is thrown.
#' 

#' @author Will Beasley
#' 
#' @examples
#' library(IalsaSynthesis) #Load the package into the current R session.
#' \dontrun{
#' path <- "./studies/eas"
#' good_name <- "u1_male_aehplus_muscle_noCog_hand_noCogSpec.out"
#' validate_filename_output(good_name, path)
#' 
#' bad_name <- "missing_something.outtttt"
#' validate_filename_output(bad_name, path)
#' }

validate_filename_output <- function( filename, path, file_extension_expected="out", underscore_count_expected=4L ) {
  testit::assert("The output filename was not specified.", 
                 !missing(filename) & !is.na(filename) & !is.null(filename))
  testit::assert("The output path was not specified.", 
                 !missing(path) & !is.na(path) & !is.null(path))
  testit::assert("The expected file extension was not specified.", 
                 !is.na(file_extension_expected) & !is.null(file_extension_expected))
  testit::assert("The expected underscore count was not specified.", 
                 !is.na(underscore_count_expected) & !is.null(underscore_count_expected))
  
  full_path <- file.path(path, filename)
  if( !file.exists(full_path) )
    stop( paste0("The output file was not found at `", full_path, "`."))
  
  extension <- tools::file_ext(full_path)
  if( extension != file_extension_expected )
    stop( paste0("The output file extension `", extension, "` did not match the expected value of `", file_extension_expected, "`."))
  
#   underscore_count_observed <- sum(grepl(pattern="_", x=filename, perl=T))
#   error_message <- paste0("An output filename should contain exactly ", underscore_count_expected, 
#     " underscores.  It had ", underscore_count_observed, ".")
#   testit::assert(error_message, underscore_count_observed==underscore_count_expected)
  
  return( invisible(TRUE) )
}
