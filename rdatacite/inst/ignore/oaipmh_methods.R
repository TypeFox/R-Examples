o_list_records <- function(url = dc_oai_base(), prefix = "oai_dc", from = NULL,
                           until = NULL, set = NULL, ...) {

  args <- dc_compact(list(verb = "ListRecords", metadataPrefix = prefix, from = from,
                          until = until, set = set))
  res <- GET(url, query = args, ...)
  stop_for_status(res)
  tt <- content(res, "text")
  xml_orig <- xml2::read_xml(tt)
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



# iterate version --------
# o_list_records <- function(url = dc_oai_base(), prefix = "oai_dc", from = NULL,
#                            until = NULL, set = NULL, ...) {
#
#   args <- dc_compact(list(verb = "ListRecords", metadataPrefix = prefix, from = from,
#                           until = until, set = set))
#   res <- GET(url, query = args, ...)
#   stop_for_status(res)
#   tt <- content(res, "text")
#   xml_orig <- xml2::read_xml(tt)
#   xml <- xml2::xml_children(xml2::xml_children(xml_orig)[[3]])
#   headers <- do.call("rbind.data.frame", dc_compact(lapply(xml, function(z) {
#     if (xml2::xml_name(z) != "resumptionToken") {
#       tmp <- xml2::xml_children(z)[[1]]
#       data.frame(rbind(unlist(xml2::as_list(tmp))), stringsAsFactors = FALSE)
#     }
#   })))
#
#   metadata <- dc_compact(lapply(xml, function(y) {
#     if (xml_name(y) != "resumptionToken") {
#       tmp <- xml2::xml_children(xml2::xml_children(xml2::xml_children(y)[[2]]))
#       lapply(tmp, function(w) {
#         as.list(setNames(xml2::xml_text(w), xml2::xml_name(w)))
#       })
#     }
#   }))
#
#   list(headers = headers, metadata = metadata)
# }
#
# iter <- 0
# token <- "characters" # define a iterator, also used for gettingn the resumptionToken
# out <- list() # define empty list to put joural titles in to
# while (is.character(token)) { # while token is class "character", keep going
#   iter <- iter + 1
#   args2 <- args
#   if (token != "characters") {
#     args2$resumptionToken <- token
#     args2$metadataPrefix <- NULL
#   }
#
#   res <- GET(url, query = args2, ...)
#   stop_for_status(res)
#   tt <- content(res, "text")
#   xml_orig <- xml2::read_xml(tt)
#   xml <- xml2::xml_children(xml2::xml_children(xml_orig)[[3]])
#   tok <- xml2::xml_text(xml2::as_list(xml[sapply(xml, xml_name) == "resumptionToken"])[[1]])
#   tok_atts <- xml2::xml_attrs(xml2::as_list(xml[sapply(xml, xml_name) == "resumptionToken"])[[1]])
#   tok <- c(token = tok, as.list(tok_atts))
#   out[[iter]] <- get_data(xml)
#   if (is(tok$token, "character")) {
#     token <- 1
#   } else {
#     token <- tok$token
#   }
# }


do.call("c", pluck(out, "headers"))

get_data <- function(x) {
  list(headers = get_headers(x), metadata = get_metadata(x))
}

get_headers <- function(x) {
  do.call("rbind_fill", dc_compact(lapply(x, function(z) {
    if (xml2::xml_name(z) != "resumptionToken") {
      tmp <- xml2::xml_children(z)[[1]]
      dat <- lapply(xml2::xml_children(tmp), function(w) {
        as.list(setNames(xml2::xml_text(w), xml2::xml_name(w)))
      })
      data.frame(rbind(unlist(dat)), stringsAsFactors = FALSE)
    }
  })))
}

get_metadata <- function(x) {
  dc_compact(lapply(x, function(y) {
    if (xml2::xml_name(y) != "resumptionToken") {
      tmp <- xml2::xml_children(y)
      status <- unlist(xml_attrs(tmp))
      if (length(status) != 0) {
        NULL
      } else {
        tmp <- xml2::xml_children(xml2::xml_children(tmp[[2]]))
        lapply(tmp, function(w) {
          as.list(setNames(xml2::xml_text(w), xml2::xml_name(w)))
        })
      }
    }
  }))
}
