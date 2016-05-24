#' @title Update the locally stored OAI-PMH data providers table.
#'
#' @description Data comes from \url{http://www.openarchives.org/Register/BrowseSites}.
#' Data includes oai-identifier (if they have one) and baes URL. The website has
#' the name of the data provider too, but not provided in the data pulled
#' down here, but you can grab the name using example below.
#'
#' @export
#' @details This table is scraped from
#' 		\url{http://www.openarchives.org/Register/BrowseSites}.
#' 		I would get it from \url{http://www.openarchives.org/Register/ListFriends},
#' 		but it does not include repository names.
#'
#' 		This function updates the table for you. Does take a while though, so
#' 		go get a coffee.
#' @param path Path to put data in.
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @seealso \code{\link{load_providers}}
#' @examples \dontrun{
#' update_providers()
#' load_providers()
#' }

update_providers <- function(path = ".", ...) {
  tt <- GET(oai_base(), ...)
  stop_for_status(tt)
  temp <- content(tt, "text", encoding = "UTF-8")
  prov <- xml2::read_html(temp)
  tab <- xml2::xml_find_all(prov, "//table")[[2]]
  children <- xml2::xml_children(tab)
  providers <- rbind.fill(lapply(children[-1], function(z) {
    data.frame(t(gsub("\n|\\s\\s+", "", xml2::xml_text(xml2::xml_children(z)[3:5]))),
               stringsAsFactors = FALSE)
  }))
  names(providers) <- c("repo_name", "base_url", "oai_identifier")
  save(providers, file = paste(path, "/", Sys.Date(), "-providers.rda", sep = ""))
}

oai_base <- function() "http://www.openarchives.org/Register/BrowseSites"
