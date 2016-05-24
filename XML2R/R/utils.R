#' Parse XML Files into XML Documents
#' 
#' Essentially a recursive call to \link{xmlParse}.
#' 
#' @param urls character vector or list of urls that point to an XML file (or anything readable by \link{xmlParse}).
#' @param async logical. Allows for asynchronous download requests. This option is passed to the \code{async} option in the \code{RCurl::getURL} function.
#' @param quiet logical. Print file name currently being parsed?
#' @importFrom plyr try_default
#' @importFrom XML xmlParse
#' @importFrom RCurl getURL
#' @importFrom RCurl url.exists
#' @export
            
urlsToDocs <- function(urls, async=TRUE, quiet=FALSE) {
  #keep only urls that exist
  urls <- urls[vapply(urls, url.exists, logical(1), USE.NAMES=FALSE)]
  if (async && length(urls) > 1) message("Performing asynchronous download. Please be patient.")
  text <- getURL(urls, async=async)
  docs <- NULL
  for (i in seq_along(text)) {
    if (!quiet) cat(urls[i], "\n")
    doc <- try_default(xmlParse(text[i], asText=TRUE), NULL, quiet = TRUE)
    if (!is.null(doc)) {
      attr(doc, "XMLsource") <- urls[i]
      docs <- c(docs, doc) #Keep non-empty documents
    }
  }
  return(docs)
}

#' Parse XML Documents into XML Nodes
#' 
#' Essentially a recursive call to \link{getNodeSet}.
#' 
#' @param docs XML documents
#' @param xpath xpath expression
#' @importFrom XML getNodeSet
#' @export

docsToNodes <- function(docs, xpath) {
  #I should really figure which class I want...
  rapply(docs, function(x) getNodeSet(x, path=xpath), 
         classes=c('XMLInternalDocument', 'XMLAbstractDocument'), how="replace")
}


#' Coerce XML Nodes into a list with both attributes and values
#' 
#' Essentially a recursive call to \link{xmlToList}.
#' 
#' @param nodes A collection of XML nodes. Should be the output from \link{docsToNodes}.
#' @importFrom XML xmlToList
#' @return A nested list with a structure that resembles the XML structure
#' @export
#' 

nodesToList <- function(nodes){
  #I should really figure which class I want...
  rapply(nodes, function(x) xmlToList(x),
    classes=c("XMLInternalElementNode", "XMLInternalNode", "XMLAbstractNode"), how="replace")
}

#' Flatten nested list into a list of observations
#' 
#' This function flattens the nested list into a list of "observations" (that is, a list of matrices with one row).
#' The names of the list that is returned reflects the XML ancestory of each observation.
#' 
#' @param l list. Should be the output from \link{nodesToList}. 
#' @param urls character vector the same length as \code{l}. Each element should map element of \code{l} to an XML file.
#' @param append.value logical. Should the XML value be appended to the observation?
#' @param as.equiv logical. Should observations from two different files (but the same ancestory) have the same name returned?
#' @param url.map logical. If TRUE, the 'url_key' column will contain a condensed url identifier (for each observation)
#' and full urls will be stored in the "url_map" element. If FALSE, the full urls are included (for each observation) 
#' as a 'url' column and no "url_map" is included.
#' @return A list where each element reflects one "observation".
#' @export

#adapted from knowledege gained from here:
#http://stackoverflow.com/questions/8139677/how-to-flatten-a-list-to-a-list-without-coercion?
##http://stackoverflow.com/questions/18862601/extract-name-hierarchy-for-each-leaf-of-a-nested-list

listsToObs <- function(l, urls, append.value=TRUE, as.equiv=TRUE, url.map=TRUE) {
  #add a prefix to the names of the list to help content to a url
  url.count <- paste0("url", seq_len(length(urls)))
  url_map <- cbind(url_key=url.count, url=urls)
  names(l) <- url.count
  #Assuming the names of the list elements holding XML values (and only values) will be NULL, we can distinguish between values/attributes
  #By also assuming that values appear immediately before their respective attributes (which appears to be the way xmlToList works),
  #we append the XML value to a row of attributes (given that their from the same node)
  name.len <- rapply(l, function(x) length(names(x)))  
  list.len <- length(name.len)
  nms <- names(name.len)
  #if (is.null(names(l[[1]]))) nms <- rep(names(l), sapply(l, length)) #overwrite names if no names (or no children) exist
  # turn '..attr' into '.attr' -- note it can show up multiple times in a name
  temp <- gsub('.attrs', 'attrs', fixed=TRUE, nms)
  idx <- gsub('.', '//', fixed=TRUE, temp)
  #select the url prefix and replace with nothing (this will give observations from different files and the same ancestory the same name)
  if (as.equiv) {
    node.sets <- sub("url([0-9]+)//", "", idx)
  } else {
    node.sets <- idx
  }
  urlz <- sub("//.*$", "", idx)
  suffix <- sub(".*//", "", node.sets)
  #tracks which XML values should be appended to the sequential row
  indicies <- which(name.len == 0 & suffix %in% "text") 
  holder <- vector('list', list.len) #placeholder for the flattened list hierarchy
  i <- 0L
  #fill up placeholder with relevant info
  rapply(l, function(x) { 
              i <<- i+1L 
              holder[[i]] <<- matrix(x, nrow=1)
              nmz <- names(x)
              if (is.null(nmz)) nmz <- "XML_value"
              colnames(holder[[i]]) <<- nmz
            })
  #Append XML_value column to the appropriate attributes
  #Note that this assumes the value always appears before the attributes in the list order
  if (append.value && length(indicies) > 0) {
    values <- holder[indicies]
    add.values <- holder[indicies+1]
    holder[indicies+1] <- mapply(function(x, y){ cbind(x, y) }, add.values, values, SIMPLIFY=FALSE)
    holder[indicies] <- NULL
    #remove the XML value elements (so they aren't duplicated)
    urlz <- urlz[-indicies]
    node.sets <- node.sets[-indicies]
  }
  if (url.map) {
    holder <- mapply(function(x, y) cbind(x, url_key=y), holder, urlz, SIMPLIFY=FALSE)
    holder[[length(holder) + 1]] <- url_map #append url_map to end of list
    node.sets <- c(node.sets, "url_map")
  } else {
    tab <- table(urlz)
    reps <- tab[match(url.count, names(tab))]
    reps[is.na(reps)] <- 0 #error handling for urls without any observations
    holder <- mapply(function(x, y) cbind(x, url=y), holder, rep(urls, reps), SIMPLIFY=FALSE)
  }
  names(holder) <- sub("//attrs", "", node.sets)
  return(holder)
}


#' Rename rows of a list
#' 
#' Sometimes, certain nodes in an XML ancestory may want to be neglected
#' before any keys are created (see \link{add_key}) or observations are aggregated (see \link{collapse}).
#' This function takes a list of "observations" (that is, a list of matrices with one row) and 
#' alters the names of that list. Note that any information lost from changing names is saved 
#' in a new column whose name is specified by \code{diff.name}.
#' 
#' @param obs list. Should be the output from \link{XML2Obs} (or \link{listsToObs}). 
#' @param namez must be equivalent to \code{names(obs)}. Intended use is to avoid unneccessarily repeating that operation.
#' @param equiv character vector with the appropriate (unique) names that should be regarded "equivalent".
#' @param diff.name character string used for naming the variable that is appended to any observations whose name was overwritten. 
#' The value for this variable is the difference in from the original name and the overwritten name.
#' @param rename.as character string to override naming of observations that are renamed.
#' @param quiet logical. Include message about how observations are being renamed?
#' @return A list of "observations". 
#' @export

re_name <- function(obs, namez, equiv, diff.name="diff_name", rename.as, quiet=FALSE){
  if (missing(equiv)) {
    warning("Must include equiv argument!")
    return(obs)
  }
  if (missing(namez)) {
    nms <- names(obs)
  } else {
    nms <- namez
  }
  if (is.null(nms)) {
    warning("The observations don't not have any names!")
    return(obs)
  }
  if (all(!equiv %in% unique(nms))) {
    warning("None of the equiv elements match the names of the observations.")
    return(obs)
  }
  baseline <- strsplit(equiv, "//")
  n <- length(equiv)
  keeps <- baseline[[1]]
  names(baseline) <- equiv
  for (i in seq_len(n-1)) { #get all the common nodes
    idx <- keeps %in% baseline[[i+1]]
    keeps <- keeps[idx]
  }
  #construct the 'label' to replace all the names that match equiv
  if (length(keeps) == 0) {
    warning("None of the XML node names seem to match. Are you sure you want to regard these names equivalent?")
    label <- "equivalent"
  } else {
    label <- paste0(keeps, collapse="//") 
  }
  if (!missing(rename.as)) label <- rename.as
  if (!quiet) message(paste0("Renaming all list elements named: \n", paste(equiv, collapse="  OR  "), "\nwith\n", label))
  diffs <- lapply(baseline, function(x) paste(x[!x %in% keeps], collapse="//")) #keeps the nodes that will be 'overwritten'
  idx <- nms %in% equiv
  #idx <- grepl(paste(equiv, collapse="||"), nms)  #grep here instead?
  #overwrite the names
  names(obs)[idx] <- label
  #get the information that was "lost" when overwritting names and append a new column accordingly
  klass <- as.character(diffs[match(nms[idx], names(diffs))]) 
  obs[idx] <- mapply(function(x, y) cbind(x, `colnames<-`(cbind(y), diff.name)), obs[idx], klass, SIMPLIFY=FALSE)
  return(obs)
}

#' Add a key to connect parents to descendants
#' 
#' This function creates a mapping from parent observations to it's descendants (which useful for merging/joining tables).
#' Either an existing value in the parent observation can be \code{recycle}d to it's descendants or a new column 
#' will be created (if \code{recycle} is missing).
#' 
#' @param obs list. Should be the output from \link{listsToObs}. 
#' @param parent character string. Should be present in the names of \code{obs}.
#' @param recycle character string that matches a variable name among \code{parent} observations.
#' @param key.name The desired column name of the newly generated key.
#' @param quiet logical. Include message about the keys being generated?
#' @return A list of observations.
#' @export

add_key <- function(obs, parent, recycle, key.name, quiet=FALSE){
  if (missing(parent)) {
    warning("You must provide the parent argument!")
    return(obs)
  }
  if (length(parent) > 1) warning("Please specify one parent at a time!")
  nms <- names(obs)
  if (is.null(nms)) {
    warning("The observations don't not have any names!")
    return(obs)
  }
  un <- unique(nms)
  if (!parent %in% un) {
    warning("The parent argument you provided does not match any observations.")
    return(obs)
  }
  fetus <- un[-which(un == parent)]
  children <- fetus[grep(paste0(parent, "//.*"), fetus)]
  if (length(children) == 0){
    warning(paste0("No children were found for the ", parent, " node."))
    return(obs)
  } else {
    if (!quiet) message(paste0("A key for the following children will be generated for the ", parent, " node:\n", paste0(children, collapse="\n")))
  }
  if (missing(recycle)) { #add an (outer) index to parent
    elders <- nms == parent
    outer_index <- seq_len(sum(elders))
    if (missing(key.name)) key.name <- "key_name"
    obs[elders] <- mapply(function(x, y) cbind(x, `colnames<-`(cbind(y), key.name)), obs[elders], outer_index, SIMPLIFY=FALSE) 
  } else { #use a specified variable for (outer) index
    elders <- nms == parent
    outer_index <- lapply(obs[elders], function(x) x[,recycle])
    if (missing(key.name)) key.name <- recycle
  }
  #now add the (inner) index for each descendent
  for (child in children) {
    kids <- nms == child
    kid.idx <- which(kids)
    relevant <- obs[elders | kids] #have to subset to relevant cases to get the proper indexing
    elder.idx <- which(names(relevant) == parent)
    timez <- diff(c(0, elder.idx))-1 #yields the number of children between each parent
    inner_index <- rep(outer_index, timez)
    obs[kid.idx] <- mapply(function(x, y) cbind(x, `colnames<-`(cbind(y), key.name)), obs[kid.idx], inner_index, SIMPLIFY=FALSE) 
  }
  return(obs)
}

#' Collapse a list of observations into a list of tables.
#' 
#' This function aggregates all observations with a similar name into a common table. Note that observations 
#' with a particular name don't need consistent variables (any missing information is filled with NAs).
#' 
#' @param obs list of observations.
#' @return Returns list with one element for each relevant XML node. Each element contains a matrix.
#' @importFrom plyr rbind.fill.matrix
#' @export

collapse_obs <- function(obs) {
  nms <- names(obs)
  #Exclude url_map from collapsing if it exists
  map <- grep("url_map", nms)
  if (length(map) != 0) {
    obs <- obs[-map]
    nms <- nms[-map] 
  }
  if (length(unique(nms)) == 1) {
    return(rbind.fill.matrix(obs))
  } else {
    return(tapply(obs, INDEX=nms, rbind.fill.matrix))
  }
}
