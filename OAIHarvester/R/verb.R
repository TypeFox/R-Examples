## OAI-PMH request issuers.

## See http://www.openarchives.org/OAI/openarchivesprotocol.html

## These can either return the result in XML form, or (by default) a
## suitable transformation thereof (see below).

oaih_get_record <-
function(baseurl, identifier, prefix = "oai_dc", transform = TRUE)
{
    request <-
        sprintf("verb=GetRecord&identifier=%s&metadataPrefix=%s",
                identifier, prefix)
    response <- OAI_PMH_issue_request(baseurl, request)
    result <- OAI_PMH_get_result(response)[["record"]]
    if(transform) oaih_transform(result) else result
        
}

oaih_identify <-
function(baseurl, transform = TRUE)
{
    response <- OAI_PMH_issue_request(baseurl, "verb=Identify")
    result <- OAI_PMH_get_result(response)
    if(transform) oaih_transform(result) else result
}

oaih_list_identifiers <-
function(baseurl, prefix = "oai_dc",
         from = NULL, until = NULL, set = NULL, transform = TRUE)
{
    request <-
        c(sprintf("verb=ListIdentifiers&metadataPrefix=%s", prefix),
          ## Note that from or until could be date/time objects ...
          if(!is.null(from)) sprintf("from=%s", from),
          if(!is.null(until)) sprintf("until=%s", until),
          if(!is.null(set)) sprintf("set=%s", set))
    OAI_PMH_gather_request_results(baseurl, request, transform)
}

oaih_list_metadata_formats <-
function(baseurl, identifier = NULL, transform = TRUE)
{
    request <-
        c("verb=ListMetadataFormats",
          if(!is.null(identifier))
          sprintf("identifier=%s", identifier))
    response <- OAI_PMH_issue_request(baseurl, request)
    result <- OAI_PMH_get_result(response)
    if(transform) oaih_transform(result) else result    
}

oaih_list_records <-
function(baseurl, prefix = "oai_dc",
         from = NULL, until = NULL, set = NULL, transform = TRUE)
{
    request <-
        c(sprintf("verb=ListRecords&metadataPrefix=%s", prefix),
          ## Note that from or until could be date/time objects ...
          if(!is.null(from)) sprintf("from=%s", from),
          if(!is.null(until)) sprintf("until=%s", until),
          if(!is.null(set)) sprintf("set=%s", set))
    OAI_PMH_gather_request_results(baseurl, request, transform)
}
         
oaih_list_sets <-
function(baseurl, transform = TRUE)
{
    request <- "verb=ListSets"
    OAI_PMH_gather_request_results(baseurl, request, transform)
}
