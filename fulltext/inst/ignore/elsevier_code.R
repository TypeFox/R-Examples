library("httr")

# Fetch -----------------------------------------------------
# doi <- "10.1016/j.ibusrev.2010.09.002" # i think not oa
doi <- "10.1016/j.tree.2013.03.008" # i think not oa
# doi <- "10.1016/j.tree.2014.05.004" # i think not oa
# doi <- "10.1016/j.celrep.2015.07.015" # oa
# url <- paste0("http://api.elsevier.com/content/article/DOI:", doi)
url <- paste0("http://api.elsevier.com/content/article/doi/", doi)
key <- Sys.getenv("CROSSREF_TDM_ELSEVIER")
# key <- getOption("crossref_tdm_elsevier")
h <- add_headers(`X-ELS-APIKey` = key, Accept = "text/xml")
# h <- add_headers(`X-ELS-APIKey` = key, Accept = "text/plain")
# p <- httr::use_proxy("136.152.208.242", 
#                      username = "scott.chamberlain", 
#                      password = "AngK#ckkLPXTHdoQKQdpjRwpNZC2k*FDe6GGEnQv")
res <- GET(url, h, verbose())
res$status_code
content(res)



# Search -----------------------------------------------------
els_search <- function(query, limit = 10, start = 1, ...) {
  args <- list(query = query, count = limit, start = start)
  hd <- add_headers(`X-ELS-APIKey` = key, Accept="application/json")
  res <- GET(els_base, query = args, hd, ...)
  stop_for_status(res)
  content(res)
}
els_base <- 'http://api.elsevier.com/content/search/index:SCIDIR'

out <- els_search(query = "cellular", config=verbose())
out$`search-results`$`opensearch:totalResults`
out$`search-results`$`opensearch:startIndex`
out$`search-results`$`opensearch:itemsPerPage`
out$`search-results`$`opensearch:Query`
out$`search-results`$link
out$`search-results`$entry
out$`search-results`$entry[[1]]
