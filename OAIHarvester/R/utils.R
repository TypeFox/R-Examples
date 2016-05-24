## Helper functions.

.get_value_of_unique_text_nodes <-
function(x, names = NULL)
{
    names(x) <- NULL
    if(is.null(names))
        names <- unique(unlist(lapply(x, xmlApply, xmlName)))
    out <- lapply(x,
                  function(n) {
                      lapply(n[names],
                             .xml_value_if_not_null,
                             NA_character_)
                  })
    matrix(unlist(out, recursive = FALSE),
           nrow = length(x), ncol = length(names),
           byrow = TRUE, dimnames = list(NULL, names))    
}

.get_value_of_any_text_nodes <-
function(x, names = NULL)
{
    if(is.null(names))
        names <- unique(unlist(lapply(x, xmlApply, xmlName)))
    ## <NOTE>
    ## We used to eventually apply
    ##    as.character(sapply(.elements_named(e, nm), xmlValue)))
    ## but apparently xmlValue() gives character() rather than an empty
    ## string for "empty" nodes ... eventually resulting in
    ##    "character(0)"
    ## entries.  Argh.  Let's use .xml_value_if_not_empty() instead.
    out <-
        do.call(cbind,
                lapply(names,
                       function(nm) {
                           lapply(x, 
                                  function(e)
                                  as.character(sapply(.elements_named(e,
                                                                      nm),
                                                      .xml_value_if_not_empty))
                                  )
                       }))
    dimnames(out) <- list(NULL, names)
    out
}

.elements_named <-
function(x, name)
    .structure(x[names(x) == name], names = NULL)

## <FIXME>
## Is the distinction between getting a single "any" node and all "any"
## nodes still useful?
## </FIXME>

.get_one_any_node <-
function(x, name)
    x[[name]][[1L]]

.get_all_any_nodes <-
function(x, name)
    lapply(.elements_named(x, name), `[[`, 1L)

## <NOTE>
## Damn.
## All of OAI-PMH is UTF-8.
## We use getURL with .encoding = "UTF-8" and the XML we get starts with
##   "<?xml version='1.0' encoding='UTF-8'?>"
## but nevertheless xmlTreeParse() and subsequent xmlValue() only return
## strings with encoding "unknown" (2010-09-15).
## The encoding argument to xmlValue() seems unused ...
## Hence, try adding the encoding back in ...

.xml_value_if_not_null <-
function(n, default)
    if(!is.null(n)) .xml_value_in_utf8(n) else default

.xml_value_if_not_empty <-
function(n)
    if(length(v <- .xml_value_in_utf8(n))) v else ""

.xml_value_in_utf8 <-
function(n)
{
    v <- xmlValue(n)
    if(is.character(v))
        Encoding(v) <- "UTF-8"
    v
}
## </NOTE>

.OAI_PMH_UTC_date_stamp <-
function(x, times_ok = TRUE)
{
    ## This could also be used in oaih_list_records(), either with a
    ## times_ok = TRUE default (users should know what they do in case
    ## they give a date/time object ...).  Alternatively, in case users
    ## give from/to and we want to be nice, we could query the
    ## repository ... 
    
    if(inherits(x, "POSIXt")) {
        ## Convert to GMT.
        x <- as.POSIXlt(x, tz = "GMT")
        x <- strftime(x, if(times_ok) "%FT%TZ" else "%F")
    }
    else {
        x <- as.character(x)
        if(!times_ok && (nchar(x) > 10L)) {
            warning("Repository only supports dates.")
            x <- substring(x, 1L, 10L)
        }
    }
    
    x
}

.structure <-
function(x, ...)
    `attributes<-`(x, c(attributes(x), list(...)))
