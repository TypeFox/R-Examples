## OAI-PMH infrastructure.

## See http://www.openarchives.org/OAI/openarchivesprotocol.html

OAI_PMH_issue_request <-
function(baseurl, request)
{
    verbose <- getOption("verbose")
    request <- paste(request, collapse = "&")

    ## <FIXME>
    ## Add support for compression eventually ...
    ## </FIXME>
    
    url <- URLencode(paste(baseurl, request, sep = "?"))

    if(verbose)
        message(gettextf("Performing request '%s'", url))

    ans <- GET()(url)
    
    ## <NOTE>
    ## E.g,
    ##   http://CRAN.R-project.org/oai
    ## redirects to
    ##   http://cran.R-project.org:8080/repo/CRANpackages
    ## One can handle such redirections using followlocation = TRUE:
    ## however, when used together with header = TRUE, this gives the
    ## headers from *all* requested pages ... giving some 3xx code.
    ## Hence, we handle the 3xx codes which indicate locations to
    ## redirect ourselves.  See e.g. "HTTP status codes 3xx" in
    ##   http://en.wikipedia.org/wiki/URL_redirection
    ## </NOTE>

    ## Look at the header first to see if we succeeded.
    h <- ans$header
    ## OAI-PMH says the Content-Type returned for OAI-PMH requests must
    ## be text/xml (even in the case of non-OK status codes?).  So let
    ## us look at the HTTP status codes directly.
    if((s <- h["Status-Code"]) != "200") {
        ## OAI-PMH says certain status codes may be useful for load
        ## balancing.
        ##   503 - Service unavailable, a Retry-After period is
        ##   specified.  Harvesters should wait this period before
        ##   attempting another OAI-PMH request.
        ##   302 - Allows the repository to temporarily redirect an
        ##   OAI-PMH request to another repository. The URI of the
        ##   temporary repository should be given by the Location field
        ##   in the HTTP response.
        ## See above for URL directions.
        if((s == "503") && !is.na(t <- h["Retry-After"])) {
            if(verbose)
                message(gettextf("Need to retry after %s seconds", t))
            Sys.sleep(as.numeric(t) + 1L)
            return(Recall(baseurl, request))
        } else if((s %in% c("300", "301", "302", "303", "307"))
                  && !is.na(l <- h["Location"])) {
             if(verbose)
                 message(gettextf("Need to redirect to %s", l))
             return(Recall(sub("[?].*", "", l), request))
        } else {
            msg <-
                sprintf("OAI-PMH request failed with HTTP status code %s",
                        s)
            txt <- h["Reason-Phrase"]
            if(!is.na(txt))
                msg <- paste(msg, sprintf("and message:\n%s", txt))
            stop(msg)
        }
    }

    ## Proceed with body.
    
    ## http://www.openarchives.org/OAI/2.0/openarchivesprotocol.htm says
    ## that the XML responses to OAI-PMH requests have the following
    ## common markup:
    ## * The first tag output is an XML declaration where the version is
    ##   always 1.0 and the encoding is always UTF-8.
    ## * The remaining content is enclosed in a root element with the
    ##   name OAI-PMH.  This element must have three attributes that
    ##   define the XML namespaces used in the remainder of the response
    ##   and the location of the validating schema.
    ## * The first two children of the root element are always
    ##   ** responseDate (a UTCdatetime indicating the time and date
    ##      that the response was sent).
    ##   ** request (indicating the protocol request that generated this
    ##      response).
    ## * The third child of the root element is either:
    ##   ** an error element that must be used in case of an error or
    ##      exception condition;
    ##   ** an element with the same name as the verb of the respective
    ##      OAI-PMH request.
    ## See http://www.openarchives.org/OAI/2.0/OAI-PMH.xsd.

    ## We will refer to the third child as the "result" in the non-error
    ## case.

    ## But what should we return?
    ## In case of an error, most likely an error condition.
    ## Hmm.  As this is the "basic" stuff, perhaps just the XML parse
    ## tree.  Let specific request issuers handle the resumptionToken
    ## accumulation as necessary ...

    nodes <- xmlTreeParse(ans$body, asText = TRUE)

    ## It would be nice to use internal nodes and xmlPathApply for
    ## more efficiently extracting nodes ... 

    result <- OAI_PMH_get_result(nodes)

    if(xmlName(result) == "error")
        stop(OAI_PMH_error(xmlAttrs(result)["code"], xmlValue(result)))

    nodes
    
}

## Get request verb from OAI-PMH request response
OAI_PMH_get_verb <-
function(nodes)
    xmlAttrs(xmlRoot(nodes)[["request"]])["verb"]

## Get result from OAI-PMH request response.
OAI_PMH_get_result <-
function(nodes)    
    xmlRoot(nodes)[[3L]]

OAI_PMH_gather_request_results <-
function(baseurl, request, transform = FALSE)
{
    ## Gather request results unless complete (no more resumption
    ## tokens) and return the aggregated "result" of the request(s).
    
    nodes <- OAI_PMH_issue_request(baseurl, request)
    ## Errors would have been thrown.
    verb <- OAI_PMH_get_verb(nodes)
    result <- OAI_PMH_get_result(nodes)
    kids <- xmlChildren(result)

    ## Even without transforming, it seems better to gather request
    ## results in a list, and combine at the end.
    chunks <- list()
    repeat {
        size <- length(kids)
        ## Assume that the resumption token comes last.
        last <- kids[[size]]
        done <- if(xmlName(last) != "resumptionToken")
            TRUE
        else {
            ## Drop resumption token from results.
            kids <- kids[-size]
            ## Done iff the resumption token is "empty".
            !length(token <- xmlValue(last))
        }
        if(transform)
            kids <- oaih_transform(kids)
        chunks <- c(chunks, list(kids))
        if(done) break
        nodes <-
            OAI_PMH_issue_request(baseurl,
                                  sprintf("verb=%s&resumptionToken=%s",
                                          verb, token))
        kids <- xmlChildren(OAI_PMH_get_result(nodes))
    }

    if(transform) {
        result <- do.call("rbind", chunks)
    } else {
        chunks <- unlist(chunks, recursive = FALSE, use.names = FALSE)
        ## Using
        ##   xmlChildren(result) <- chunks
        ## adds names and hence duplicates storage ...
        result$children <- chunks
    }

    result    
}

OAI_PMH_error <-
function(code, info, call = NULL)
{
    msg <- sprintf("Received condition '%s'", code)
    if(length(info) && nzchar(info))
        msg <- paste(msg, sprintf("with diagnostic:\n%s", info))
    err <- list(code = code, info = info, message = msg, call = call)
    class(err) <- c("OAI_PMH_error", "error", "condition")
    err
}
