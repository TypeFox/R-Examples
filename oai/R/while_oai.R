while_oai <- function(url, args, token, as, dumper=NULL, dumper_args=NULL, ...) {
  iter <- 0
  token <- "characters"
  out <- list()
  while (is.character(token)) {
    iter <- iter + 1
    args2 <- args
    if (token != "characters") {
      args2$resumptionToken <- token
      args2$from <- NULL
      args2$until <- NULL
      args2$set <- NULL
      args2$metadataPrefix <- NULL
    }

    res <- GET(url, query = args2, ...)
    stop_for_status(res)
    tt <- content(res, "text", encoding = "UTF-8")

    # try parsing
    parsed <- try( read_xml_safely(tt), silent=TRUE )
    if(inherits(parsed, "try-error")) {
      parsed <- try( xml2::read_html(tt), silent=TRUE )
      warning("read_xml parsing failed, but read_html succeeded")
      if( inherits(parsed, "try-error") ) {
        fname <- tempfile(
          pattern=paste0("oaidump_", chartr(" :", "_-", Sys.time())),
          tmpdir=".", fileext="xml")
        cat(tt, file=fname)
        stop(paste0("cannot parse downloaded XML, dumped raw XML to file ", fname))
      }
      is_html <- TRUE
    } else {
      is_html <- FALSE
    }

    handle_errors(parsed)
    tok <- get_token(parsed, verb=args2$verb, is_html=is_html)
    # `as` determines what the `dumper` gets
    if (as == "raw") {
      res <- tt
    } else {
      if(is_html) {
        warning("malformed XML - keeping raw text even though `as` is ", dQuote(as))
        res <- tt
      } else {
        xml_verb <- xml2::xml_children(xml2::xml_children(parsed)[[3]])   # TODO xpath here
        res <- switch(args$verb,
               ListRecords = get_data(xml_verb, as = as),
               ListIdentifiers = parse_listid(xml_verb, as = as),
               ListSets = get_sets(xml_verb, as = as)
        )
      }
    }
    # Collect values returned by `dumper` if they are not NULL
    if (is.null(dumper)) {
      out[[iter]] <- res
    } else {
      valid_dumper(dumper, dumper_args)
      dumper_res <- do.call("dumper", c(list(res = res, args = args2, as = as), dumper_args))
      if (!is.null(dumper_res))
        out[[iter]] <- dumper_res
    }
    if (tok$token == "") {
      token <- 1
    } else {
      token <- tok$token
    }
  }

  switch(args$verb,
         ListRecords = do.call("c", out),
         ListIdentifiers = do.call("c", out),
         ListSets = out
  )
}

# Get resumptionToken from parsed XML/HTML
# @param x result of read_xml or read_html
# @param is_html is it XML or HTML (incl. malformed XML)
get_token <- function(x, verb, is_html=FALSE) {
  xp <- paste0("/*[local-name()='OAI-PMH']/*[local-name()='", verb, "']/*[local-name()='",
               if(is_html) "resumptiontoken" else "resumptionToken", "']" )
  node <- xml2::xml_find_all(x, xp)
  if(length(node) == 0) {
    return( list(token="") )
  }
  else {
    if(length(node) > 1) {
      warning("more than one match - using last")
      node <- node[length(node)]
    }
    return( c(
      token=xml2::xml_text(node),
      as.list(xml2::xml_attrs(node))
    ) )
  }
}
