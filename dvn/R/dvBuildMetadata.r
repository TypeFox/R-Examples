dvBuildMetadata <- function(..., format='dcterms'){
    if(!format=='dcterms')
        stop('Only format=="dcterms" is currently supported')
    
    pairls <- list(...)
    dublincore <- c(
        "abstract","accessRights","accrualMethod","accrualPeriodicity",
        "accrualPolicy","alternative","audience","available",
        "bibliographicCitation","conformsTo","contributor","coverage",
        "created","creator","date","dateAccepted","dateCopyrighted",
        "dateSubmitted","description","educationLevel","extent","format",
        "hasFormat","hasPart","hasVersion","identifier","instructionalMethod",
        "isFormatOf","isPartOf","isReferencedBy","isReplacedBy","isRequiredBy",
        "issued","isVersionOf","language","license","mediator","medium",
        "modified","provenance","publisher","references","relation","replaces",
        "requires","rights","rightsHolder","source","spatial","subject",
        "tableOfContents","temporal","title","type","valid")
    if(any(!names(pairls) %in% dublincore))
        stop('All names of `...` parameters must be in Dublin Core')
    if(!'title' %in% names(pairls))
        stop('"title" is a required metadata field')
    entry <- newXMLNode('entry', namespaceDefinitions=
                c(  'http://www.w3.org/2005/Atom',
                    dcterms='http://purl.org/dc/terms/'))
    dcchild <- function(x,y)
        dcnode <- newXMLNode(y, x, namespace='dcterms')
    addChildren(entry, mapply(dcchild,pairls,names(pairls)))
    entry <- paste('<?xml version="1.0" encoding="UTF-8" ?>\n',toString.XMLNode(entry),sep='')
    class(entry) <- c(class(entry),'dvMetadata')
    attr(entry,'formatName') <- format
    return(entry)
}
