#' 	Perform Confidence Run
#'
#' 	This function starts the `confidence tool`. The results will be stored 
#'     in a subdirectory in the current working directory. See details section
#'     below.
#'     
#'
#' 	@param x name of the input file or a \code{\link{data.frame}} containing the 
#'     input. If \code{x = \link{NULL}} (the default) a file dialog will 
#'     appear for interactive selection of an input file. See the package 
#'     vignette for details about the file format. 
#'	@param tmpdir directory to store temporary files (for debugging only)
#'  @param browse load resulting report directly in a browser? \code{TRUE} or \code{FALSE}
#'  
#'  @details This function will create a subdirectory
#'      \itemize{
#'          \item{in the same directory as \code{x}, in case \code{x} is a filename or}
#'          \item{in the current working directory (see \code{\link{getwd}}), in case 
#'     \code{x} is a \code{\link{data.frame}}}.
#'    }
#'    The computer should have write permission to this directory, if not
#'    an error message will be raised. The
#'    subdirectory contains an HTML-report with all analysis results. 
#'    For convenience, the results are also stored in CSV-format (tables) 
#'    and png-format (figures) for further processing.
#'    
#'  @seealso \link{confidence} and the package vignette 
#'      (\code{vignette("confidence")}).
#'
#' 	@export
conf <-
function(x = NULL, tmpdir = tempfile(pattern = "confidence"), browse = TRUE) {

    # prevent potential problems with dates in other locales
    oldLocale <- Sys.getlocale("LC_TIME")
    on.exit(Sys.setlocale("LC_TIME", oldLocale))
    Sys.setlocale("LC_TIME", "C")

    # read data
    if (!is.data.frame(x)) {

        # interactive selection of filename
        if (is.null(x)) {
            x <- tk_choose.files(
                default = "", 
                caption = "Select input file",
                multi = FALSE, 
                filters = matrix(data = c("confidence input", ".csv"), nrow = 1)
            )
        }
    
        # check if filename exists
        if (length(x) == 0L || !file.exists(x)) {
            stop(
                sprintf("Input file %s does not exist", sQuote(x)), 
                call. = FALSE
            )
        }
        
        # store valid filename
        filename <- x

        # read file data
    	x <- try(read.csv(file = filename, as.is = TRUE), silent = TRUE)
    	if (inherits(x, "try-error")) {
    		stop(
                "Errors occurred while reading input file ", sQuote(filename),
    			call. = FALSE
    		)
    	}
    }

    # initialization message
    message("The confidence tool is running...")

    # convert data to instance of class 'conf_input'
    x <- conf_input(x)

    # create output directory
    if (exists("filename")) {
        outputDir <- dirname(filename)
    } else {
        outputDir <- getwd()
    }
    isWritable <- file.access(names = outputDir, mode = 2) == 0L
    if (!isWritable) {
        stop(sprintf("no write permission for directory:\n'%s'", 
                     outputDir), call. = FALSE)
    }
    outputDir <- file.path(
        normalizePath(outputDir), 
        paste0("output", format(Sys.time(), format = "%Y%m%dT%H%M%S"))
    )
    
    # create HTML-file
    write_html(x = x, outputDir = outputDir, browse = browse)
    
    # store csv file
    write.csv(
        x = as.data.frame(mya(x)), 
        file = file.path(outputDir, "output.csv"), 
        row.names = FALSE
    )

    # finalization
    message("The confidence run has been completed successfully.")
 }


#'  @inheritParams conf
#'  
#'  @rdname conf
#'  
#'  @export
CONF <-
function(x = NULL, tmpdir = tempfile(pattern = "confidence"), browse = TRUE) {
    conf(x = x, tmpdir = tmpdir, browse = browse)
}