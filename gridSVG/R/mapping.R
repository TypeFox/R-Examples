testUniqueMappings <- function(x) {
    idNodes <- getNodeSet(x, "//*[@id]")
    ids <- sapply(idNodes, function(x) xmlGetAttr(x, "id"))
    length(ids) == length(unique(ids))
}

formatTypeMapping <- function(x, type) {
    objs <- x[x$type == type, c("name", "suffix", "selector", "xpath")]
    if (! nrow(objs))
        return(NULL)
    objNames <- unique(objs$name)
    objList <- vector("list", length(objNames))
    names(objList) <- objNames
    for (i in 1:length(objNames))
        objList[[i]] <- objs[objs$name == objNames[i], c("suffix", "selector", "xpath")]
    objList
}

formatMappings <- function(x) {
    list(vps = formatTypeMapping(x, "vp"),
         grobs = formatTypeMapping(x, "grob"),
         refs = formatTypeMapping(x, "ref"),
         id.sep = getSVGoption("id.sep"),
         prefix = get("prefix", envir = .gridSVGEnv))
}

exportMappings <- function(x) {
    x <- formatMappings(x)
    paste("var gridSVGMappings = ", toJSON(x), ";\n", sep = "")
}

gridSVGMappingsGen <- function() {
    mappings <- NULL
    function(newmappings = NULL) {
        if (is.null(newmappings)) {
            mappings
        } else {
            mappings <<- newmappings
        }
    }
}

gridSVGMappings <- gridSVGMappingsGen()

getSVGMappings <- function(name, type, result = "id") {
    if (! type %in% c("vp", "grob", "ref"))
        stop("'type' must be one of 'vp', 'grob' or 'ref'")
    if (! result %in% c("id", "selector", "xpath"))
        stop("'result' must be one of 'id', 'selector' or 'xpath'")
    # Because the list itself uses vp/grob, rewrite
    type <- paste(type, "s", sep = "")
    mappings <- gridSVGMappings()
    if (is.null(mappings))
        stop("gridSVGMappings() must be initialised")
    nameData <- mappings[[type]][[name]]
    if (is.null(nameData))
        stop("Name not found")

    if (result == "id")
        paste0(mappings$prefix, name, mappings$id.sep, nameData$suffix)
    else
        nameData[[result]]
}

readMappingsJS <- function(filename) {
  jsData <- readLines(filename)
  jsData <- gsub("var gridSVGMappings = ", "", jsData)
  jsonData <- gsub(";$", "", jsData)
  fromJSON(paste0(jsonData, collapse = "\n"))
}
