fieldInfo <- rbind(
    c("pathData",           "factor",    "additional"),
    c("pathDescription",    "factor",    "additional"),
    c("identification",     "character", "obligatory"),
    c("fullName",           "character", "optional"),
    c("dirName",            "character", "obligatory"),
    c("files",              "character", "obligatory"),
    c("originalColsNumber", "integer",   "obligatory"),
    c("originalColsNames",  "character", "optional"),
    c("originalColsType",   "character", "obligatory"),
    c("delete",             "character", "optional"),
    c("responsePos",        "integer",   "obligatory"),
    c("attributesNumber",   "integer",   "computed"),
    c("attributesType",     "character", "computed"),
    c("responseType",       "character", "computed"),
    c("cases",              "integer",   "obligatory"),
    c("available",          "logical",   "additional"))

readDSListFromXML <- function(filename)
{
    dsList <- xmlToDataFrame(filename, homogeneous=FALSE)

    for (i in which(fieldInfo[, 3] == "optional")) {
        if (is.null(dsList[[fieldInfo[i, 1]]])) {
            dsList[[fieldInfo[i, 1]]] <- ""
        }
        dsList[is.na(dsList[[fieldInfo[i, 1]]]), fieldInfo[i, 1]] <- ""
    }

    for (i in which(fieldInfo[, 2] == "integer" & fieldInfo[, 3] == "obligatory")) {
        dsList[[fieldInfo[i, 1]]] <- as.integer(dsList[[fieldInfo[i, 1]]])
    }

    attributesNumber <- rep(NA, times=nrow(dsList))
    attributesType <- rep(NA, times=nrow(dsList))
    responseType <- rep(NA, times=nrow(dsList))

    for (i in seq.int(length.out=nrow(dsList))) {
        originalColsType <- split.comma(dsList$originalColsType[i])
        delete <- as.integer(split.comma(dsList$delete[i]))
        attributesNumber[i] <- dsList$originalColsNumber[i] - length(delete) - 1L
        attributesType[i] <- paste(originalColsType[ - c(delete, dsList$responsePos[i])], collapse=",")
        responseType[i] <- originalColsType[dsList$responsePos[i]]
    }
    dsList$attributesNumber <- attributesNumber
    dsList$attributesType <- attributesType
    dsList$responseType <- responseType

    dsList <- dsList[, fieldInfo[fieldInfo[, 3] != "additional", 1]]
    dsList
}

saveDSListAsXML <- function(dsList, filename)
{
    stopifnot(is.data.frame(dsList))
    cnames <- names(dsList)
    datasets <- list()
    for (i in 1:nrow(dsList))
    {
        components <- list()
        for (j in which(fieldInfo[, 3] %in% c("obligatory", "optional"))) {
            if (dsList[i, j] != "")
                components[[length(components)+1]] <- xmlNode(cnames[j], xmlTextNode(dsList[i,j]))
        }
        datasets[[length(datasets)+1]] <- xmlNode("dataset", .children=components)
    }
    X <- xmlNode("contents", .children=datasets)
    saveXML(X, file=filename)
}

