get_technical_information <-
function(x)
{
    ## Helper function for getting a single BibTeX entry.
    get_single_BibTeX_entry <- function(x) {
        ## First, get the type.
        type <- .jcall(x, "Lweka/core/TechnicalInformation$Type;",
                       "getType")
        type <- .jcall(type, "S", "toString")
        ## Next, get the fields.
        fields <- .jcall(x, "Ljava/util/Enumeration;", "fields")
        ## <FIXME>
        ## For some reason,
        ##   fields <- x$fields()
        ## does not work.
        ## </FIXME>
        ## And try looping ...
        tags <- values <- character()
        while(.jcall(fields, "Z", "hasMoreElements")) {
            field <- .jcast(.jcall(fields,
                                   "Ljava/lang/Object;",
                                   "nextElement"),
                            "weka/core/TechnicalInformation$Field")
            ## Now, field is a weka.core.TechnicalInformation$Field
            ## object.
            ## Get its display string.
            tags <- c(tags, .jcall(field, "S", "toString"))
            ## And its value.
            values <- c(values, .jcall(x, "S", "getValue", field))
        }
        names(values) <- tags
        BibTeX_entry(type, .jcall(x, "S", "getID"), values)
    }

    ## Now get all the BibTeX entries.
    ## First, instantiate and get the technical information.
    info <- .jcall(x, "Lweka/core/TechnicalInformation;",
                   "getTechnicalInformation")
    ## The first entry.
    entries <- list(get_single_BibTeX_entry(info))
    ## And possibly additional ones.
    additional <- .jcall(info, "Ljava/util/Enumeration;", "additional")
    while(.jcall(additional, "Z", "hasMoreElements")) {
        a <- .jcall(additional, "Ljava/lang/Object;", "nextElement")
        entries <- c(entries, list(get_single_BibTeX_entry(a)))
    }

    BibTeX_db(list = entries)
}

BibTeX_entry <-
function(type, key = "", fields)
    .structure(list(type = type, key = key, fields = fields),
               class = "BibTeX_entry")

BibTeX_db <-
function(..., list = NULL)
{
    ## For the time being, assume that we are given lists of BibTeX
    ## entries.
    `class<-`(c(list(...), list), "BibTeX_db")
}

format.BibTeX_entry <-
function(x, offset = 0L, ...)
{
    prefix <- paste(rep.int(" ", offset), collapse = "")
    fields <- x$fields
    c(sprintf("%s@%s{%s,", prefix, toupper(x$type), x$key),
      strwrap(sprintf("%s = {%s},", names(fields), fields),
              indent = offset + 2L, exdent = offset + 4L),
      sprintf("%s}", prefix))
}

format.BibTeX_db <-
function(x, offset = 0L, ...)
{
    n <- length(x)
    if(!n) return(character())
    unlist(mapply(c, rep.int(list(""), n),
                  lapply(x, format.BibTeX_entry, offset = offset)))[-1L]
}

print.BibTeX_db <-
function(x, ...)
{
    writeLines(format(x))
    invisible(x)
}

print.BibTeX_entry <-
function(x, ...)
{
    writeLines(format(x))
    invisible(x)
}
