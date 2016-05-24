.checkPlaFrame <- function(data = NULL,
                           design = "",
                           latin = FALSE,
                           blocks = FALSE,
                           crd = FALSE)
{
    OK <- TRUE
    design <- .string2design(design)
    if (is.null(data)) {
        warning("No data!")
        OK <- FALSE
    }
    if (max(dim(data)) < 1) {
        warning("No data!")
        OK <- FALSE
    }
    Names <- dimnames(data)[[2]]
    if (!any(Names == "Response")) {
        warning("'Response' not found in data")
        OK <- FALSE
    }
    if (!any(Names == "Dilution")) {
        warning("'Response' not found in data")
        OK <- FALSE
    }
    if (!any(Names == "Sample")) {
        warning("'Sample' not found in data")
        OK <- FALSE
    }
    if (blocks | (design == "rbd"))
        if (!any(Names == "Replicate")) {
            warning(paste0("'Replicate' not found in data for design ", design, "[rbd]"))
            OK <- FALSE
        }
    if (latin | (design == "lsd"))
        if (!any(Names == "Row")) {
            warning(paste0("'Row' not found in data for design ", design, "[lsd]"))
            OK <- FALSE
        }
    if (latin | (design == "lsd"))
        if (!any(Names == "Column")) {
            warning(paste0("'Column' not found in data for design ", design, "[lsd]"))
            OK <- FALSE
        }
    return(OK)
}
