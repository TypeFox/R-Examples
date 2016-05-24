oaih_size <-
function(baseurl, from = NULL, until = NULL, set = NULL)
{
    request <-
        c("verb=ListIdentifiers&metadataPrefix=oai_dc",
          ## Note that from or until could be date/time objects ...
          if(!is.null(from)) sprintf("from=%s", from),
          if(!is.null(until)) sprintf("until=%s", until),
          if(!is.null(set)) sprintf("set=%s", set))

    nodes <- tryCatch(OAI_PMH_issue_request(baseurl, request),
                      error = identity)
    if(inherits(nodes, "error")) {
        if(grepl("noRecordsMatch", conditionMessage(nodes)))
            return(0)
        else
            stop(nodes)
    }
    
    result <- OAI_PMH_get_result(nodes)
    kids <- xmlChildren(result)
    size <- length(kids)
    last <- kids[[size]]
    if(xmlName(last) == "resumptionToken") {
        size <- xmlAttrs(last)["completeListSize"]
        if(is.null(size)) size <- NA_real_
    }

    as.numeric(size)
}
