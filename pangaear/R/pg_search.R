#' Search the Pangaea database
#'
#' @export
#' @param query (character) Query terms. You can refine a search by prefixing the term(s) with a
#' category, one of citation, reference, parameter, event, project, campaign, or basis.
#' See examples.
#' @param count (integer) Number of items to return.
#' @param env (character) Type of data to search, one of "all", "sediment", "water", "ice", "atomosphere"
#' @param bbox  (numeric) A bounding box, of the form: minlon, minlat, maxlon, maxlat
#' @param mindate,maxdate (character) Dates to search for, of the form "2014-10-28"
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @return data.frame
#' @details This is a thin wrapper around the GUI search interface on the page
#' \url{https://www.pangaea.de}. Everything you can do there, you can do here.
#' @examples \dontrun{
#' pg_search(query='water')
#' pg_search(query='water', count=2)
#' pg_search(query='water', count=20)
#' pg_search(query='water', env="water")
#' pg_search(query='water', env="sediment")
#' pg_search(query='water', mindate="2013-06-01", maxdate="2013-07-01")
#' pg_search(query='water', bbox=c(-124.2, 41.8, -116.8, 46.1))
#' pg_search(query='citation:Archer')
#' pg_search(query='reference:Archer')
#' pg_search(query='parameter:"carbon dioxide"')
#' pg_search(query='event:M2-track')
#' pg_search(query='event:TT011_2-CTD31')
#' pg_search(query='project:Joint Global Ocean Flux Study')
#' pg_search(query='campaign:M2')
#' pg_search(query='basis:Meteor')
#' }

pg_search <- function(query, count=10, env="all", bbox=NULL, mindate=NULL, maxdate=NULL, ...) {
  check_if(count, c("numeric", "integer"))
  check_if(env, "character")
  check_if(mindate, "character")
  check_if(maxdate, "character")
  args <- pgc(list(count = count, q = query, env = capwords(env), mindate = mindate, maxdate = maxdate))
  if (!is.null(bbox)) args <- c(args, as.list(setNames(bbox, c('minlon', 'minlat', 'maxlon', 'maxlat'))))
  res <- GET(sbase(), query = args, ...)
  stop_for_status(res)
  html <- read_html(content(res, "text", encoding = "UTF-8"))
  nodes <- xml_find_all(html, "//li")
  dat <- lapply(nodes, parse_res)
  as_data_frame(do.call("rbind.data.frame", lapply(dat, as_data_frame)))
}

parse_res <- function(x){
  doi <- sub("https?://doi.pangaea.de/", "", xml_attr(xml_find_all(x, './/p[@class="citation"]/a'), "href"))
  citation <- xml_text(xml_find_all(x, './/p[@class="citation"]/a'))
  tab <- xml_find_all(x, './/table/tr')
  #supp <- xml_text(xml_find_all(tab[[1]], ".//td")[2])
  supp <- xml_text(xml_find_one(xml_parent(xml_find_all(tab, ".//td[contains(.,'Supplement')]")), './/td[@class="content"]'))
  size <- strextract(xml_text(xml_find_all(xml_parent(xml_find_all(tab, ".//td[contains(.,'Size')]")), './/td[@class="content"]')), "[[:digit:]]+")
  #score <- strextract(strsplit(xml_text(xml_find_all(tab[[3]], ".//td")), "Score:")[[1]][2], "[[:digit:]]+\\.[[:digit:]]+")
  score <- strextract(strsplit(xml_text(xml_find_all(tab, './/td[@class="datasetid"]')), "Score:")[[1]][2], "[[:digit:]]+\\.[[:digit:]]+")
  lis <- list(doi = doi, score = as.numeric(score), size_datasets = as.numeric(size),
       citation = citation, supplement_to = supp)
  lis[vapply(lis, length, 1) == 0] <- NA
  lis
}

check_if <- function(x, cls) {
  if (!is.null(x)) {
    if (!class(x) %in% cls) {
      stop(substitute(x), " must be of class: ", paste0(cls, collapse = ", "), call. = FALSE)
    }
  }
}
