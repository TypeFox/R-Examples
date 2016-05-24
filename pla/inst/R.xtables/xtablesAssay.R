xtablesAssay <- function(data = NULL,
                         model = NULL,
                         fits = NULL,
                         frame = NULL,
                         imputeMissing = FALSE,
                         potency.digits = 4,
                         showRelative = TRUE, ...) {
    selectFun <- function (array) NULL
    if(!is.null(model))
        selectFun <- model@selectFun
    if(!is.null(data)) {
        if (class(data) == "assayFrame") {
            if(is.null(frame)) {
                xtableHead(data, ...)
                if (.string2design(Data@design) == "lsd")
                    xtableLatinSquare(data@tableRaw)
            }
        } else {
            xtableTable(data, sampleLabels = fits@inpArgs$sampleLabels,
                        selectFun = selectFun, ...)
        }
    }
    if(!is.null(frame)) {
        if(!is.null(data))
            xtableHead(data, ...)
        design <- "crd"
        if (!is.null(model))
            design <- .string2design(model@design)
        if (design == "lsd") {
            ## Works only for latin squares !!!
            xtableLatinSquare(frame)
        } else
            xtableTableCS(frame, data)
    }
    if (imputeMissing) {
        ## if(!is.null(data))
        ##     xtableHead(data, ...)
        if (model@design == "lsd") {
            ## Works only for latin squares !!!
            xtableLatinSquare(model@data)
        } else
            xtableTableCS(model@data, data, ...)
    }
    if(!is.null(model)) {
        cat(paste0("\\medskip \n \\noindent \n"))
        xtableModel(model)
    }
    xtableAnova(fits)
    xtableValidity(fits)
    xtableParameters(fits)
    xtablePotency(fits, digits = potency.digits, showRelative = showRelative)
}
