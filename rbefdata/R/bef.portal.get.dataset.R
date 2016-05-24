#' Fetch primary data from a BEFdata portal dataset.
#'
#' This function fetches data associated with a BEFdata portal dataset. By
#' default it will fetch the CSV file of a dataset. You need to provide the
#' function with a dataset id which you can find in the URL of the dataset on
#' the BEFdata portal. As this usually requires authentication you need to set
#' your credentials in the options (bef.options("user_credentials" = "asdfpoj").
#' You can find the credentials inside of your profile page on the BEFdata portal.
#' The credentials ensure you have the rights to download the data.
#'
#' The function returns a dataset object which you can store in a variable as
#' shown in the examples below. The object also offers additional information
#' by attributes. You can query the information via the attributes() function
#' which is also shown in the examples. If you like to fetch multiple datasets
#' you can use the apply functions provided by R see example below.
#'
#' @param id This is the ID of a dataset on a BEFdata portal.
#' @param curl If the function is used inside a loop, call getCurlHandle() first
#'        and pass in the returned value here. This avoids an unnecessary footprint.
#' @param \dots Arguments passed to \code{\link[RCurl]{getURLContent}}
#'
#' @return The function returns a data frame of the dataset. An error is thrown when the dataset is
#'         not found or if you don't have the rights to access it.
#'
#' @examples \dontrun{
#'         datset1 = bef.portal.get.dataset(id=8)
#'         metadata1 = attributes(dataset1)
#'
#'         curl = getCurlHandle()
#'         ids = c(8,70)
#'         dataset_list = lapply(ids, function(x) bef.portal.get.dataset_by(id = x, curl = curl))
#'         metadata = attributes(dataset_list[[1]])
#'       }
#' @import RCurl
#' @export bef.portal.get.dataset bef.get.dataset bef.get.dataset_by bef.portal.get.dataset_by
#' @aliases bef.get.dataset bef.get.dataset_by bef.portal.get.dataset_by

bef.portal.get.dataset <-  bef.get.dataset <- bef.get.dataset_by <- bef.portal.get.dataset_by <- function(id, curl=getCurlHandle(), ...) {
  dataset_url = dataset_url(id, user_credentials= bef.options("user_credentials"), type = "csv2")
  response_body = getURLContent(dataset_url, curl = curl, ...)
  if(getCurlInfo(curl)$response.code != 200) {
    msg = sprintf("Dataset(id=%d) not found or not accessible. Please check your credentials and make sure you have access right for it.", id)
    stop(msg)
  }
  dataset = read.csv(text = response_body)
  metadata = bef.portal.get.metadata(id)
  attributes(dataset) = c(attributes(dataset), metadata)
  return(dataset)
}

#' Fetch a list of datasets for a keyword
#'
#' This function fetches a list of datasets associated with a BEFdata portal keyword.
#'
#' @param keyword The keyword you like to fetch the associated datasets for
#'
#' @examples \dontrun{
#'         list = bef.portal.get.datasets.for_keyword(keyword = "carbon")
#'       }
#' @import RCurl
#' @import rjson
#' @export bef.portal.get.datasets.for_keyword bef.get.datasets.for_keyword
#' @aliases bef.get.datasets.for_keyword

bef.portal.get.datasets.for_keyword <- bef.get.datasets.for_keyword <- function(keyword) {
  keyword_json = fromJSON(getURL(paste0(bef.options('url'),"/keywords.json")))
  names = unlist(lapply(keyword_json, function(x) (x$name)))
  ids = unlist(lapply(keyword_json, function(x) (x$id)))
  keyword_summary = data.frame(key = names, id = ids)

  if(!missing(keyword)) {
    grep_matches <- c(keyword)
    matches <- unique(grep(paste(grep_matches,collapse="|"), keyword_summary$key, value=TRUE))
    position = which(keyword_summary$key %in% matches)
    get_keyword_id = keyword_summary$id[position]
    keyword_datasets_api_list = (unique(unlist(lapply(get_keyword_id, function(x) paste0(keyword_url(x),".csv")))))
    dataset_list = keyword_datasets_api_list
    dataset_info = lapply(dataset_list, function(x) read.csv(x)[,1:2])
    titles = unique(unlist(lapply(dataset_info, function(x) x$title)))
    ids = unique(unlist(lapply(dataset_info, function(x) x$id)))
    id_title_df = data.frame(id = ids, title = titles)
    if (dim(id_title_df)[1] == 0 ) {
      print("Sorry no datasets are tagged with this keyword!")
    } else {
      return(id_title_df)
    }
  }
}
