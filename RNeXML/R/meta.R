## Utilities for adding additional metadata


#' Constructor function for metadata nodes
#' 
#' @param property specify the ontological definition together with it's namespace, e.g. dc:title
#' @param content content of the metadata field
#' @param rel Ontological definition of the reference provided in href 
#' @param href A link to some reference
#' @param datatype optional RDFa field
#' @param id optional id element (otherwise id will be automatically generated).  
#' @param type optional xsi:type.  If not given, will use either "LiteralMeta" or "ResourceMeta" as 
#'  determined by the presence of either a property or a href value. 
#' @param children Optional element containing any valid XML block (XMLInternalElementNode class, see the XML package for details).  
#' @details User must either provide property+content or rel+href.  Mixing these will result in potential garbage. 
#' The datatype attribute will be detected automatically from the class of the content argument.  Maps from R class
#' to schema datatypes are as follows: 
#' character - xs:string, 
#' Date - xs:date,
#' integer - xs:integer,
#' numeric - xs:decimal,
#' logical - xs:boolean
#' 
#' @examples
#' meta(content="example", property="dc:title")
#' @export 
#' @seealso \code{\link{nexml_write}}
#' @include classes.R
meta <- function(property = character(0), 
                 content = character(0), 
                 rel = character(0), 
                 href = character(0), 
                 datatype = character(0), 
                 id = character(0),
                 type = character(0),
                 children = list()){
  if(is.logical(content))
    datatype <- "xsd:boolean"
  else if(is(content, "Date"))
    datatype <- "xsd:date"
  else if(is.numeric(content))
    datatype <- "xsd:decimal"
  else if(is.character(content))
    datatype <- "xsd:string"
  else if(is.integer(content))
    datatype <- "xsd:integer"
  else 
    datatype <- "xsd:string"

  # Having assigned the datatype, 
  # the content text must be written as a string
  content <- as.character(content)


  if(length(id) == 0)
    id <- nexml_id("m")
   
  if(is(children, "XMLAbstractNode") || is(children, "XMLInternalNode"))
    children <- list(children)

  if(length(property) > 0){ ## avoid 
    if(is.null(content) && length(children) == 0) ## Avoid writing when content is missing, e.g. prism:endingpage is blank
      NULL
    else
      new("meta", content = content, datatype = datatype, 
          property = property, id = id, 'xsi:type' = "LiteralMeta",
          children = children)
  } else if(length(rel) > 0){
    if(is.null(href))
      NULL
    else
      new("meta", rel = rel, href = href, 
          id = id, 'xsi:type' = "ResourceMeta",
          children = children)
  } else {
    new("meta", content = content, datatype = datatype, 
        rel = rel, href = href, id = id, 'xsi:type' = type,
        children = children)
  }
}


## Common helper functions 


nexml_citation <- function(obj){
  if(is(obj, "BibEntry"))
    class(obj) <- "bibentry"
  if(is(obj, "bibentry")){
    out <- lapply(obj, function(obj){
      if(length(grep("--", obj$pages)) > 0){
        pgs <- strsplit(obj$pages, "--")[[1]]
        start_page <- pgs[[1]]
        end_page <- if(length(pgs)>1) pgs[[2]] else " "
      } else if(length(grep("-", obj$pages)) > 0){
        pgs <- strsplit(obj$pages, "-")[[1]]
        start_page <- pgs[[1]]
        end_page <- if(length(pgs)>1) pgs[[2]] else " "
      } else {
        start_page <- NULL
        end_page <- NULL
      }
      list_of_metadata_nodes <- plyr::compact(c(list(
        meta(content=obj$volume, 
            property="prism:volume"),
        meta(content=obj$journal, 
            property="dc:publisher"),
        meta(content=obj$journal, 
            property="prism:publicationName"),
        meta(content = end_page, 
            property="prism:endingPage"),
        meta(content=start_page, 
            property="prism:startingPage"),
        meta(content=obj$year, 
            property="prism:publicationDate"),
        meta(content=obj$title,
            property="dc:title")),
        lapply(obj$author, function(x){
        meta(content = format(x, c("given", "family")),
             property="dc:contributor") 
        })))
        citation_elements = new("ListOfmeta", list_of_metadata_nodes)
        meta(content=format(obj, "text"), 
            property="dcterms:bibliographicCitation",
            children = lapply(citation_elements, as, "XMLInternalElementNode"))
    })
    out 
  }
}



#' Concatenate meta elements into a ListOfmeta
#' 
#' Concatenate meta elements into a ListOfmeta
#' @param x,... meta elements to be concatenated, e.g. see \code{\link{meta}}
#' @param recursive  logical, if 'recursive=TRUE', the function 
#' descends through lists and combines their elements into a vector.
#' @return a listOfmeta object containing multiple meta elements. 
#' @examples 
#' c(meta(content="example", property="dc:title"),
#'   meta(content="Carl", property="dc:creator"))
#' 
setMethod("c", 
          signature("meta"),
          function(x, ..., recursive = FALSE){
            elements <- list(x, ...)
#            if(recursive)
            elements <- meta_recursion(elements)
            new("ListOfmeta", elements)

          })


#' Concatenate ListOfmeta elements into a ListOfmeta
#' 
#' Concatenate ListOfmeta elements into a ListOfmeta
#' @param x,... meta or ListOfmeta elements to be concatenated, e.g. see \code{\link{meta}}
#' @param recursive  logical, if 'recursive=TRUE', the function 
#' descends through lists and combines their elements into a vector.
#' @return a listOfmeta object containing multiple meta elements. 
#' @include classes.R
#' @examples 
#' metalist <- c(meta(content="example", property="dc:title"),
#'               meta(content="Carl", property="dc:creator"))
#' out <- c(metalist, metalist) 
#' out <- c(metalist, meta(content="a", property="b")) 
setMethod("c", 
          signature("ListOfmeta"),
          function(x, ..., recursive = FALSE){
            elements <- list(x, unlist(...))
            elements <- meta_recursion(elements)
            new("ListOfmeta", elements)

          })




meta_recursion <- function(elements){
  i <- 1
  out <- vector("list")
  for(e in elements){
    if(length(e) > 0){
      if(is(e, "meta")){
        out[[i]] <- e
        i <- i + 1
      } else if(is.list(e)){
        out <- c(out, meta_recursion(e))
        i <- length(out) + 1 
      }
    }
  }
out
}
