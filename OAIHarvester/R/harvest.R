oaih_harvest <-
function(baseurl,
         prefix = "oai_dc", from = NULL, until = NULL, set = NULL,
         transform = TRUE)
{
    ## Harvest the OAI repository at baseurl.

    id <- oaih_identify(baseurl)
    ## Provides info on compression and granularity.
    ## Currently, compression is unused, but if 'gzip' is possible, we
    ## should be able to use gzcon(url(.....)) in the OAI_PMH request.

    set <- if(!is.null(set)) {
        ## Find out if the repository provides sets.
        sets <- tryCatch(oaih_list_sets(baseurl), error = identity)
        ## Match the given sets against the available sets.
        if(inherits(sets, "error"))
            stop("Repository does not provide sets.")
        sets <- unlist(sets[, "setSpec"])
        bad <- ! set %in% sets
        if(any(bad))
            stop(sprintf(ngettext(length(bad),
                                  "Set %s is unavaible.",
                                  "Sets %s are unavailable."),
                         paste(sQuote(set[bad]), collapse = ", ")))
        as.list(set)
    } else {
        list(NULL)
    }

    ## If 'from' or 'until' are given, check against the repository's
    ## granularity.
    times_ok <- id$granularity == "YYYY-MM-DDThh:mm:ssZ"
    if(!is.null(from))
        from <- .OAI_PMH_UTC_date_stamp(from, times_ok)
    if(!is.null(until))
        until <- .OAI_PMH_UTC_date_stamp(until, times_ok)

    if(!identical(prefix, "oai_dc")) {
        ## Determine available metadata formats.
        formats <-
            unlist(oaih_list_metadata_formats(baseurl)[,
                                                       "metadataPrefix"])
        ## If 'prefix' is given, match against the available formats.
        if(!is.null(prefix)) {
            bad <- ! prefix %in% formats
            if(any(bad))
                stop(sprintf(ngettext(length(bad),
                                      "Metadata prefix %s is unavaible.",
                                      "Metadata prefixes %s are unavailable."),
                             paste(sQuote(prefix[bad]), collapse = ", ")))
        } else {
            prefix <- formats
        }
    }

    ## Now get records:
    do.call(rbind,
            Map(function(p, s)
                tryCatch(oaih_list_records(baseurl, p, from, until, s,
                                           transform),
                         error = function(e) NULL),
                rep.int(prefix, length(set)),
                rep.int(set, rep.int(length(prefix), length(set)))))
}
