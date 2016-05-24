library("httr")
library("xml2")

o_list_records <- function(url = dc_oai_base(), prefix = "oai_dc", from = NULL,
                           until = NULL, set = NULL, ...) {

  args <- dc_compact(list(verb = verb, metadataPrefix = prefix, from = from,
                          until = until, set = set))
  res <- GET(url, query = args, ...)
  stop_for_status(res)
  tt <- content(res, "text")
  xml_orig <- xml2:::read_xml(tt)
  xml <- xml2::xml_children(xml2::xml_children(xml_orig)[[3]])
  headers <- do.call("rbind.data.frame", dc_compact(lapply(xml, function(z) {
    if (xml2::xml_name(z) != "resumptionToken") {
      tmp <- xml2::xml_children(z)[[1]]
      data.frame(rbind(unlist(xml2::as_list(tmp))), stringsAsFactors = FALSE)
    }
  })))

  metadata <- dc_compact(lapply(xml, function(y) {
    if (xml_name(y) != "resumptionToken") {
      tmp <- xml2::xml_children(xml2::xml_children(xml2::xml_children(y)[[2]]))
      lapply(tmp, function(w) {
        as.list(setNames(xml2::xml_text(w), xml2::xml_name(w)))
      })
    }
  }))

  list(headers = headers, metadata = metadata)
}

# args <- list(verb = "ListRecords", metadataPrefix = 'oai_dc',
#              from = '2011-06-01T', until = '2011-07-01T')
# res <- GET(dc_oai_base(), query = args)
# tt <- content(res, "text")
# xml_orig <- xml2:::read_xml(tt)
# xml <- xml2::xml_children(xml2::xml_children(xml_orig)[[3]])
# headers <- do.call("rbind.data.frame", dc_compact(lapply(xml, function(z) {
#   if (xml2::xml_name(z) != "resumptionToken") {
#     tmp <- xml2::xml_children(z)[[1]]
#     data.frame(rbind(unlist(xml2::as_list(tmp))), stringsAsFactors = FALSE)
#   }
# })))
#
# metadata <- dc_compact(lapply(xml, function(y) {
#   if (xml_name(y) != "resumptionToken") {
#     tmp <- xml2::xml_children(xml2::xml_children(xml2::xml_children(y)[[2]]))
#     lapply(tmp, function(w) {
#       as.list(setNames(xml2::xml_text(w), xml2::xml_name(w)))
#     })
#   }
# }))

# list(headers = headers, metadata = metadata)

# xml3 <- XML::xmlParse(tt)
# getNodeSet(xml3, "/ListRecords")
# request <- c(sprintf("verb=ListRecords&metadataPrefix=%s",
#                      prefix), if (!is.null(from)) sprintf("from=%s", from),
#              if (!is.null(until)) sprintf("until=%s", until), if (!is.null(set)) sprintf("set=%s",
#                                                                                          set)

# xml <- xmlParse(tt)
# xmlRoot(xml)[[3L]]
# xml
# getNodeSet(xml, "//ListSets")

# nodes <- OAIHarvester:::OAI_PMH_issue_request(baseurl, "verb=ListSets")
# verb <- OAIHarvester:::OAI_PMH_get_verb(nodes)
# kids <- xmlChildren(OAIHarvester:::OAI_PMH_get_result(nodes))


# chunks <- list()
# repeat {
#   size <- length(kids)
#   last <- kids[[size]]
#   done <-
#   if (xmlName(last) != "resumptionToken") {
#     TRUE } else {
#     kids <- kids[-size]
#     !length(token <- xmlValue(last))
#   }
#   if (transform)
#     kids <- OAIHarvester:::oaih_transform(kids)
#   chunks <- c(chunks, list(kids))
#   if (done)
#     break
#   nodes <- OAIHarvester:::OAI_PMH_issue_request(baseurl, sprintf("verb=%s&resumptionToken=%s",
#                                                   verb, token))
#   kids <- xmlChildren(OAIHarvester:::OAI_PMH_get_result(nodes))
# }
# if (transform) {
#   result <- do.call("rbind", chunks)
# }
# else {
#   chunks <- unlist(chunks, recursive = FALSE, use.names = FALSE)
#   result$children <- chunks
# }
# result
# #
#
# baseurl <- "http://epub.wu.ac.at/cgi/oai2"; oaih_list_sets(baseurl)
