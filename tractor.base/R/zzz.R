#' @import methods ore reportr
#' @importFrom grDevices col2rgb colorRamp dev.cur dev.off gray heat.colors rainbow
#' @importFrom graphics image layout lines locator par plot strwidth text
#' @importFrom stats na.omit
NULL

.onLoad <- function (libname, pkgname)
{
    if (is.null(getOption("tractorFileType")))
    {
        fileType <- toupper(Sys.getenv("TRACTOR_FILETYPE"))
        if (isTRUE(fileType %in% .FileTypes$typeNames))
            options(tractorFileType=as.vector(fileType))
        else
            options(tractorFileType="NIFTI_GZ")
    }
    
    if (is.null(getOption("tractorOutputPrecision")))
    {
        outputPrecision <- tolower(Sys.getenv("TRACTOR_OUTPUT_PRECISION"))
        if (isTRUE(outputPrecision %in% c("single","double")))
            options(tractorOutputPrecision=as.vector(outputPrecision))
        else
            options(tractorOutputPrecision="double")
    }
}
