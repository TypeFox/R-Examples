#  -----------------------------------------------------------------------------
#  pyOptions
#  =========
#' @title Options for the PythonInR package
#' @description a function for getting and setting options in the PythonInR package.
#' @param option a character giving the option to set or get. 
#' @param value the new value of the option.
#' @details The following options are available:
#' \itemize{
#'   \item \strong{numpyAlias} a character giving the numpy alias, the 
#'         default value is "numpy".
#'   \item \strong{useNumpy} a logical giving if numpy should be used
#'         if getting and setting matrices.
#'   \item \strong{pandasAlias}  a character giving the pandas alias, the 
#'         default value is "pandas".
#'   \item \strong{usePandas} a logical giving if pandas should be used
#'         if getting and setting data.frames.
#'   \item \strong{winPython364} a logical indicating if Python 3 64-bit
#'         under windows is used, this option is set automatically at startup
#'         and shouldn't be changed.
#' }
#' @examples
#' \dontrun{
#' pyOptions()
#' pyExec("import numpy as np")
#' pyOptions("numpyAlias", "np")
#' pyOptions("useNumpy", TRUE)
#' pyExec("import pandas as pd")
#' pyOptions("pandasAlias", "pd")
#' pyOptions("usePandas", TRUE)
#' }
#  -----------------------------------------------------------------------------   
pyOptions <-
local({
    options <- list(numpyAlias="numpy", useNumpy=FALSE, pandasAlias="pandas",
                    usePandas=FALSE, winPython364=FALSE)
    function(option, value) {
        if (missing(option)) return(options)
        if (missing(value)){
            options[[option]]
        }else{
            if ( class(value) != class(options[[option]]) ){
                stop(sprintf('"%s" has to be "%s"', option, class(options[[option]])))
            }
            options[[option]] <<- value
        }
    }
})

