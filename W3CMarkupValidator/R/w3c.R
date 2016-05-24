w3c_markup_validate_baseurl_default <-
function()
    Sys.getenv("W3C_MARKUP_VALIDATOR_BASEURL",
               "http://validator.w3.org/check?")

w3c_markup_validate_baseurl <-
local({
    baseurl <- w3c_markup_validate_baseurl_default()
    function(new) {
        if(missing(new))
            baseurl
        else {
            baseurl <<- if(is.null(new))
                w3c_markup_validate_baseurl_default()
            else
                new
        }
    }
})

w3c_markup_validate_message_elements <-
    c("line", "col", "message", "messageid", "explanation", "source")

w3c_markup_validate <-
function(baseurl = w3c_markup_validate_baseurl(),
         uri = NULL, file = NULL, string = NULL,
         opts = list())
{
    out <- list(valid = TRUE,
                errorcount = 0L, errors = NULL,
                warningcount = 0L, warnings = NULL)
    class(out) <- "w3c_markup_validate"

    ## Be nice and add the question mark at the end if not given.
    if(substring(baseurl, nchar(baseurl)) != "?")
        baseurl <- paste0(baseurl, "?")

    h <- new_handle()
    handle_setopt(h, .list = opts)
    
    response <- if(!is.null(uri)) {
        curl_fetch_memory(sprintf("%suri=%s;output=soap12", baseurl,
                                  URLencode(uri)),
                          h)
    } else {
        if(!is.null(file)) {
            handle_setform(h,
                           uploaded_file = form_file(file, "text/html"),
                           output = "soap12")
        } else if(!is.null(string)) {
            handle_setform(h,
                           fragment = paste(string, collapse = "\n"),
                           output = "soap12")
        } else
            stop("You must specify one of 'uri', 'file' or 'string'.")
        curl_fetch_memory(baseurl, h)
    }

    ## See <http://validator.w3.org/docs/api.html>.

    status <- response$status_code
    fields <- .status_and_headers(response)

    if((as.integer(status) %/% 100) != 2L)
        stop(sprintf("Validation request failed with status %s and message:\n%s",
                     status, trimws(fields["Reason-Phrase"])))

    status <- fields["X-W3C-Validator-Status"]
    if(is.na(status))
        stop("Not a W3C validator.")
    valid <- status == "Valid"
    if(!valid && (status != "Invalid"))
        stop("Validation failed.")
    ## Maybe provide status code or message in this case?
    out$valid <- as.logical(valid)

    ## Note that Valid is equivalent to no errors.
    errorcount <- as.integer(fields["X-W3C-Validator-Errors"])
    warningcount <- as.integer(fields["X-W3C-Validator-Warnings"])

    ## If there are no errors or warnings, there is no point looking at
    ## the body (as long as we do not obtain debug info).
    if((errorcount == 0L) && (warningcount == 0L)) return(out)

    doc <- read_xml(rawToChar(response$content))
    ns <- xml_ns(doc)
    nodes <- xml_find_all(doc, "/env:Envelope/env:Body", ns)
    if(length(nodes) != 1L)
        stop("Result format does not appear to be SOAP.")
    doc <- nodes[[1L]]

    messages <- function(doc, path) {
        elements <- w3c_markup_validate_message_elements
        x <- lapply(xml_find_all(doc, path, ns),
                    function(e)
                    trimws(.xml_children_values(e))[elements])
        ## In case we found no messages ...
        if(!length(x)) return(NULL)
        y <- as.data.frame(do.call(rbind, x), stringsAsFactors = FALSE)
        colnames(y) <- elements
        y        
    }

    if(errorcount > 0L) {
        out$errorcount <- errorcount
        out$errors <-
            messages(doc, "//m:errors/m:errorlist/m:error")
    }

    if(warningcount > 0L) {
        out$warningcount <- warningcount
        out$warnings <-
            messages(doc, "//m:warnings/m:warninglist/m:warning")
    }

    out
}

print.w3c_markup_validate <-
function(x, ...)
{
    writeLines(format(x, details = FALSE))
    invisible(x)
}

format.w3c_markup_validate <-
function(x, details = TRUE, ...)
{
    fmt <- function(m) {
        sprintf("  %s %s  %s",
                format(c("line", m[, "line"]), justify = "right"),
                format(c(" col", m[,  "col"]), justify = "right"),
                c("message",
                  gsub("[[:space:]]*\n[[:space:]]*", "  ",
                       m[, "message"])))
    }

    y <- sprintf("Valid: %s (errors: %d, warnings: %d)",
                 x$valid,
                 x$errorcount,
                 x$warningcount)

    details <- if(is.character(details)) {
        !is.na(pmatch(c("e", "w"), details))
    } else
        rep_len(details, 2L)
    
    if(details[1L] && (NROW(x$errors) > 0L)) {
        y <- c(y, "Errors:", fmt(x$errors))
    }
    if(details[2L] && (NROW(x$warnings) > 0L)) {
        y <- c(y, "Warnings:", fmt(x$warnings))
    }

    y
}

as.data.frame.w3c_markup_validate <-
function(x, row.names = NULL, optional = FALSE, ...)
{
    y <- matrix("", 0L, length(w3c_markup_validate_message_elements) + 1L)
    y <- rbind(as.data.frame(y, stringsAsFactors = FALSE),
               if(NROW(m <- x$errors)) {
                   cbind(m, category = "error")
               },
               if(NROW(m <- x$warnings)) {
                   cbind(m, category = "warning")
               })
    names(y) <- c(w3c_markup_validate_message_elements, "category")
    y
}

w3c_markup_validate_db <-
function(x, names = NULL)
{
    if(inherits(x, "w3c_markup_validate")) {
        x <- list(x)
    }
    if(is.null(names)) {
        names <- names(x)
        if(is.null(names))
            names <- as.character(seq_along(x))
    }
    names(x) <- names
    ## Better put on the names first: otherwise handling failures is
    ## more tricky.
    ind <- sapply(x, inherits, "error")
    if(any(ind)) {
        failures <- x[ind]
        x <- x[!ind]
        attr(x, "failures") <- failures
    }
    if(!all(sapply(x, inherits, "w3c_markup_validate")))
        stop("all elements must be w3c_markup_validate results")
    class(x) <- "w3c_markup_validate_db"
    x
}

`[.w3c_markup_validate_db` <-
function(x, i)
{
    y <- NextMethod("[")
    class(y) <- class(x)
    y
}

print.w3c_markup_validate_db <-
function(x, ...)
{
    writeLines(sprintf("Valid: %d out of %d (errors: %d, warnings: %d)",
                       sum(vapply(x, `[[`, NA, "valid"), na.rm = TRUE),
                       length(x),
                       sum(vapply(x, `[[`, 0L, "errorcount")),
                       sum(vapply(x, `[[`, 0L, "warningcount"))))
    if(n <- length(attr(x, "failures"))) {
        writeLines(sprintf("Failures: %d", n))
    }
    invisible(x)
}

as.data.frame.w3c_markup_validate_db <-
function(x, row.names = NULL, optional = FALSE, ...)
{
    y <- lapply(x, as.data.frame)
    lens <- sapply(y, NROW)
    do.call(rbind,
            Map(cbind,
                name = Map(rep.int, names(x), lens),
                y,
                stringsAsFactors = FALSE))
}

c.w3c_markup_validate_db <-
function(..., recursive = FALSE)
{
    args <- list(...)
    ## Could check for all elements to inherit from
    ## w3c_markup_validate_db.
    ## Note that the following does not necessarily give unique names:
    ## we could ensure these be prefixing with (possibly generated)
    ## argument names.
    w3c_markup_validate_db(unlist(args, recursive = FALSE),
                           unlist(lapply(args, names)))
}

w3c_markup_validate_files <-
function(files, baseurl = w3c_markup_validate_baseurl(), opts = list()) {
    verbose <- interactive()
    results <-
        lapply(files,
               function(f) {
                   if(verbose) message(sprintf("Processing %s ...", f))
                   tryCatch(w3c_markup_validate(baseurl = baseurl,
                                                file = f,
                                                opts = opts),
                            error = identity)
               })
    w3c_markup_validate_db(results, files)
}

w3c_markup_validate_uris <-
function(uris, baseurl = w3c_markup_validate_baseurl(), opts = list()) {
    verbose <- interactive()
    results <-
        lapply(uris,
               function(u) {
                   if(verbose) message(sprintf("Processing %s ...", u))
                   tryCatch(w3c_markup_validate(baseurl = baseurl,
                                                uri = u,
                                                opts = opts),
                            error = identity)
               })
    w3c_markup_validate_db(results, uris)
}

inspect <-
function(x, ...)
    UseMethod("inspect")

inspect.w3c_markup_validate <-
function(x, details = TRUE, ...)
{
    writeLines(format(x, details = details))
}

inspect.w3c_markup_validate_db <-
function(x, details = TRUE, full = FALSE, ...)
{
    if(!full) {
        counts <-
            sapply(x, function(e) e$errorcount + e$warningcount)
        x <- x[counts > 0L]
    }
    nmx <- names(x)
    writeLines(paste(unlist(Map(function(h, b)
                                paste(c(paste0(h, ":"),
                                        format(b, details = details)),
                                      collapse = "\n  "),
                                nmx, x)),
                     collapse = "\n\n"))
}

is_invalid <-
function(x)
    identical(x$valid, FALSE)

.status_and_headers <-
function(response)
{
    lines <- parse_headers(response$headers)
    parts <- strsplit(lines[1L], " ")[[1L]]
    s <- c("Status-Code" = parts[2L],
           "Reason-Phrase" = paste(parts[-c(1L, 2L)], collapse = " "))
    lines <- lines[-1L]
    h <- sub("[^:]+: (.*)", "\\1", lines)
    names(h) <- sub("([^:]+):.*", "\\1", lines)
    c(s, h)
}

.xml_children_values <-
function(x)
{
    kids <- xml_children(x)
    y <- xml_text(kids)
    names(y) <- xml_name(kids)
    y
}
