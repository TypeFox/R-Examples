#' @import XML stringr rjson
NULL

## Data ###

#' Cached version of daten.berlin.de dataset feed.
#' Stored as XML tree in text form
#' 
#' @docType data
#' @keywords datasets
#' @source \url{http://daten.berlin.de/datensaetze/rss.xml}
#' @format XML tree with timestamp attribute 'last_updated'
#' @name cached_datasets_feed
NULL

# globalVariables call to appease R CMD check
globalVariables('cached_datasets_feed')

## Generic functions ##

#' Gets metadata for a dataset or list of datasets
#' @param where where to look for metadata
#' @param ... optional additional arguments to methods
#' @export
#' @examples
#' # query
#' \dontrun{
#' result <- searchBerlinDatasets(query = "vornamen")
#' summary(result)
#' 
#' # options to get metadata:
#' # pass in URL
#' metadata <- getDatasetMetaData(result[[2]]$link)
#' # returns 'berlin_data_dataset' object with list of available resources
#' class(metadata)
#' 
#' # pass in object
#' metadata <- getDatasetMetaData(result[[2]])
#' # same result as passing in link
#' 
#' # pass in berlin_data_list with query results to get metadata 
#' # for all objects in list
#' metadata <- getDatasetMetaData(result)
#' summary(metadata)
#' summary(metadata[[1]])
#' }
getDatasetMetaData <- function(where, ...) {
  UseMethod("getDatasetMetaData")
}

#' Downloads a resource
#' @param x a resource, or object containing resources, which can be downloaded 
#' @param ... optional additional arguments to methods, e.g.:
#'   \itemize{
#'     \item{message.on.fail}{logical: show message on download failure?}
#'     \item{message.on.succeed}{logical: show message on download success?}
#'   }
#' @export
#' @examples
#' 
#' # query, select a dataset, get its metadata
#' \dontrun{
#' result <- searchBerlinDatasets(query = "vornamen")
#' summary(result)
#' dataset <- getDatasetMetaData(result[[2]]) 
#' summary(dataset)
#' 
#' # pick a resource to download
#' resource <- resources(dataset)[[23]]
#' data <- download(resource)
#' # returns a data.frame
#' class(data)
#' # or download multiple resources simultaneously
#' data <- download(resources(dataset)[1:6])
#' # returns a list of data.frames
#' class(data)
#' class(data[[1]])
#' 
#' # download all resources in dataset
#' data <- download(dataset)
#' 
#' # turn off individual notifications for failed downloads 
#' # due to unsupported file formats, URL schemes, etc.
#' data <- download(dataset, message.on.fail=FALSE)
#' 
#' # turn off all individual notifications for downloads
#' data <- download(dataset, message.on.fail=FALSE, message.on.succeed=FALSE)
#' 
#' # pass in other arguments to download function
#' result <- searchBerlinDatasets(query = "stolpersteine")
#' dataset <- getDatasetMetaData(result[[2]])
#' resource_list <- resources(dataset)
#' data <- download(resource_list[[1]])
#' # gives wrong output, so we try a different argument for 'sep'
#' data <- download(resource_list[[1]], sep=',')
#' }
download <- function(x, ...) {
  UseMethod("download")
}

#' Gets the resources from an object
#' @param object an object with resources
#' @param ... optional additional arguments to methods
#' @export
resources <- function(object, ...) {
  UseMethod("resources")
}

## Functions for main usage ##

#' Queries daten.berlin.de
#' 
#' 
#' @param query a query string to search daten.berlin.de
#' @param ... optional additional arguments
#' @export
#' @examples
#' \dontrun{
#' result <- searchBerlinDatasets(query = "stolpersteine")
#' summary(result)
#' dataset <- getDatasetMetaData(result[[2]])
#' summary(dataset)
#' resource_list <- resources(dataset)
#' summary(resource_list)
#' data <- download(resource_list[[1]], sep=',')
#' }
#' 
searchBerlinDatasets <- function(query, ...) {
  stopifnot(length(query) == 1 && is.character(query))
  result <- searchData(query, ...) 
  result
}

#' Parses and downloads the meta data for a dataset
#' @param dataset_url the url of the dataset
#' @export
parseMetaData <- function(dataset_url) {
  stopifnot(is.character(dataset_url))
  stopifnot(length(dataset_url) == 1)
  parsed_data <- htmlParse(dataset_url)
  title_nodeset <- xpathApply(parsed_data,  "//h1[@id='page-title']")
  title <- str_trim(xmlValue(title_nodeset[[1]]))
  resources_nodeset <- getNodeSet(parsed_data, "//div[@class='dataset_ressource']")
  stripTags <- function(htmlString) {
    text <- gsub("<.*>", "", htmlString)
    text <- gsub("<.*>", "", htmlString)
    text <- gsub("\\n", "", text)
    str_trim(text)
  }
  resources_list <- lapply(resources_nodeset, function(res) {
    sub_doc <- xmlDoc(res)
    title_doc <- getNodeSet(sub_doc, "//div[@class='resource_notes']/div")
    if (length(title_doc) == 1) {
      title <- stripTags(xmlValue(title_doc[[1]]))
    } else {
      title <- "<could not parse title>"
    }
    if (nchar(title) == 0) {
      title <- "<no title>"
    }
    field_labels <- getNodeSet(sub_doc, "//div[@class='field-label']")
    field_items <- getNodeSet(sub_doc, "//div[@class='field-items']")
    cleaned_field_labels <- unlist(lapply(field_labels, function(l)str_trim(xmlValue(l))))
    cleaned_field_items <- lapply(field_items, function(l)str_trim(xmlValue(l)))
    findIndex <- function(search_key)grep(search_key, cleaned_field_labels, ignore.case = TRUE)
    data_format_index <- findIndex("format")
    hash_index <- findIndex("hash")
    url_index <- findIndex("url")
    language_index <- findIndex("sprache")
    indexExists <- function(index)length(index) == 1 && index > 0
    if (indexExists(url_index) && indexExists(hash_index) 
        && indexExists(data_format_index)
        && indexExists(language_index)) {
      res_url <- unlist(xpathApply(xmlDoc(field_items[[url_index]]),  
                                   "//a[@href]", xmlGetAttr, "href"))
      structure(list(
        title = title,
        url = res_url,
        hash = cleaned_field_items[[hash_index]],
        format = cleaned_field_items[[data_format_index]],
        language = cleaned_field_items[[language_index]]
      ), class = "berlin_data_resource")
    } else {
      NULL
    }
  })
  
  structure(list(
    title = title,
    resources = structure(Filter(function(r)!is.null(r), resources_list), 
                          class = "berlin_data_resource_list")
  ), class = "berlin_data_dataset")
}

# Searches through daten.berlin
#
# param query the query string
# param xml.url the url to the rss feed
# usage Internal
searchData <- function(query, 
                       xml.url = "http://daten.berlin.de/datensaetze/rss.xml") {
  datasets_feed <- tryCatch(xmlParse(xml.url),
                           error=function(e) {
                             if (xml.url=="http://daten.berlin.de/datensaetze/rss.xml") {
                               message('Could not establish connection to daten.berlin.de')
                               data(cached_datasets_feed, envir = environment())
                               message('Using stored list of Berlin datasets')
                               message(paste0('Last updated: ', attr(cached_datasets_feed, 'last_updated')))
                               datasets_feed = xmlParse(cached_datasets_feed)
                             } else {
                               stop(e)
                             }
                           })
  items <- getNodeSet(datasets_feed, "//item")
  cleaned_items <- lapply(items, function(item) {
    structure(
      list(
        title = xmlValue(getNodeSet(item, "title")[[1]]),
        description = xmlValue(getNodeSet(item, "description")[[1]]),
        link = xmlValue(getNodeSet(item, "link")[[1]]),
        pub_date = xmlValue(getNodeSet(item, "pubDate")[[1]]),
        creator = gsub('^.*">|</a>$', '',
                       xmlValue(getNodeSet(item, "dc:creator", namespaces=c(dc="http://purl.org/dc/elements/1.1/"))[[1]]))
      ),
      class = "berlin_data_dataset_info"
    )
  })
  filtered_items <- Filter(function(item) {
    search_result <- grep(query, c(item$title, item$description), 
                          ignore.case = TRUE)
    length(search_result) > 1 || search_result > 0
  }, cleaned_items)
  if (length(filtered_items) > 0) {
    structure(
      filtered_items,
      class = "berlin_data_list"
    )
  } else {
    structure(list(query = query), class = "berlin_data_query_no_results")
  }
}

#' Updates cached_datasets_feed
#' 
#' @export
#' @param xml.url the url to the rss feed
updateCachedFeed <- function(xml.url="http://daten.berlin.de/datensaetze/rss.xml") {
  datasets_feed <- xmlParse(xml.url)
  cached_datasets_feed <- saveXML(datasets_feed)
  attr(cached_datasets_feed, 'last_updated') <- date()
  save(cached_datasets_feed, file='data/cached_datasets_feed.rda')
  message(paste0('Updated data/cached_datasets_feed.rda. Last version: ', attr(cached_datasets_feed, 'last_updated')))
}
