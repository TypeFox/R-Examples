xtabularAssay <- function(data = NULL,
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
                        tabularx = TRUE,
                        selectFun = selectFun, ...)
        }
        xtableLabels(data, select = "2*")
    }
    if(!is.null(frame)) {
        if(!is.null(data))
            xtableHead(data, tabularx = TRUE, ...)
        design <- "crd"
        if (!is.null(model))
            design <- .string2design(model@design)
        if (design == "lsd") {
            ## Works only for latin squares !!!
            xtableLatinSquare(frame)
        } else
            xtableTableCS(frame)
        xtableLabels(data, select = "3*")
    }
    if (imputeMissing) {
        ## if(!is.null(data))
        ##     xtableHead(data, tabularx = TRUE, ...)
        if (model@design == "lsd") {
            ## Works only for latin squares !!!
            xtableLatinSquare(model@data)
        } else
            xtableTableCS(model@data, data, ...)
        xtableLabels(data, select = "4*")
    }
    if(!is.null(model)) {
        cat(paste0("\\medskip \n \\noindent \n"))
        xtableModel(model, tabularx = TRUE)
        xtableLabels(data, select = "5*")
    }
    xtableAnova(fits)
    xtableLabels(data, select = "6*")
    xtableValidity(fits)
    xtableLabels(data, select = "7*")
    xtableParameters(fits)
    xtableLabels(data, select = "8*")
    xtablePotency(fits, digits = potency.digits, showRelative = showRelative)
    xtableLabels(data, select = "9*")
    xtableLabels(data, select = "10*")
}
