## parse docsums (esummary) ####
parse_esummary <- function(object) {
  retmode <- object$retmode()
  if (!object$no_errors()) {
    warning("Errors parsing DocumentSummary", call. = FALSE)
    return(list()) 
  } else if (retmode == 'json') {
    return(jsonlite::fromJSON(object$get_content("json")))
  } else if (retmode == 'xml') {
    x <- object$get_content("xml")
    nodes <- xset(x, '/eSummaryResult/DocSum')
    if (length(nodes) != 0) {
      uids <- xvalue(x, '/eSummaryResult/DocSum/Id')
    } else {
      nodes <- xset(x, '/eSummaryResult/DocumentSummarySet/DocumentSummary')
      uids <- vapply(nodes, XML::xmlGetAttr, name = "uid", FUN.VALUE = "")
    }
  }
  docsum <- {
    docsum_list <- lapply(nodes, parse_docsum) 
    flattened_docsum <- flatten2(docsum_list)
    # check if all docsums have same number of tags
    if (length(unique(vapply(flattened_docsum, length, 0))) > 1L) {
      warning("DocSum records have a different numbers of tags.", call. = FALSE)
      setNames(flattened_docsum, uids)
    } else {
      data.frame(Id = uids, do.call("rbind", flattened_docsum), stringsAsFactors = FALSE)
    }
  }
  docsum
}

# Parse a DocSum recursively and return it as a named list
parse_docsum <- function(ds) {
  parsefun <- docsum_parser(XML::xmlName(ds))
  parsefun(ds)
}

docsum_parser <- function(version) {
  switch(version,
         DocSum = function(ds) {
           items <- XML::xmlChildren(ds)
           items <- items[names(items) == "Item"]
           value <- vector("list", length(items))
           for (i in seq_along(items)) {
             isize <- unname(XML::xmlSApply(items[[i]], XML::xmlSize))
             if (length(isize) == 0L || all(isize == 0L)) {
               value[[i]] <- trimws(XML::xmlValue(items[[i]])) %|char|% NA
             } else {
               value[[i]] <- Recall(items[[i]])
             }
           }
           names(value) <- vapply(items, XML::xmlGetAttr, name = "Name", FUN.VALUE = "", USE.NAMES = FALSE)
           value
         },
         DocumentSummary = function(ds) {
           items <- XML::xmlChildren(ds)
           value <- vector("list", length(items))
           for (i in seq_along(items)) {
             isize <- unname(XML::xmlSApply(items[[i]], XML::xmlSize))
             if (length(isize) == 0L || all(isize == 0L)) {
               value[[i]] <- trimws(XML::xmlValue(items[[i]])) %|char|% NA
             } else {
               value[[i]] <- Recall(items[[i]])
             }
           }
           names(value) <- lapply(items, XML::xmlName)
           value
         })
}

