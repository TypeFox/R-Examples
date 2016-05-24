readSDML <- function(file = "", text = NULL, validate = FALSE,
                     read.description = FALSE, ...)
{
    txt <- file == "" && !is.null(text)
    tree <- xmlTreeParse(file = if(txt) text else file,
                         asText = txt,
                         validate=validate,
                         handlers=handlersSDML(), asTree=TRUE,
                         ignoreBlanks = FALSE, trim = FALSE, ...)

    ## is this file a StatDataML file ?
    if (xmlName(xmlRoot(tree)) != "StatDataML")
        stop("This is no StatDataML file");

    ## get root node
    statxml <- xmlRoot(tree)

    ## read dataset
    dataset <- readDatasetSDML(statxml[["dataset"]])

    ## if requested, handle description
    if (read.description) {
        ## read it
        description <- readDescriptionSDML(statxml[["description"]])

        if(!is.null(description)) {
            ## if there is a description only, use an empty list as a dummy object to carry it
            if(is.null(dataset)) dataset <- list()

            ## set description as an attribute
            attr(dataset, "SDMLdescription") <- description
        }
    }

    return(dataset)
}
