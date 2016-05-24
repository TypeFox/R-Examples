parse_IETF_language_tag <-
function(x, expand = FALSE)
{
    n <- length(x)
    y <- rep.int(list(character()), n)
    names(y) <- x

    ## How nice should we be?
    ## Allow for empty or missing elements ...

    pos <- seq_along(x)
    if(any(ind <- (is.na(x) | (x == "")))) {
        pos <- pos[!ind]
        x <- x[pos]
    }

    ## See <http://www.ietf.org/rfc/rfc4646.txt>.
    ## Language tags can be of the form (in ABNF, see
    ## <http://tools.ietf.org/rfc/rfc4234.txt>): 
    ##   langtag / privateuse / grandfathered
    ## where
    ##   privateuse    = ("x"/"X") 1*("-" (1*8alphanum))
    ##   grandfathered = 1*3ALPHA 1*2("-" (2*8alphanum))

    re_privateuse <- "[xX]((-[[:alnum:]]{1,8}){1,})"

    ## Grandfathered tags must really be determined by exact matching.

    ind <- !is.na(match(x,
                        IANA_language_subtag_registry_grandfathered_table$Tag))
    if(any(ind)) {
        y[pos[ind]] <- as.list(sprintf("Grandfathered=%s", x[ind]))
        x[ind] <- ""
        pos <- pos[!ind]
    }

    if(length(pos)) {    
        pat <- sprintf("^%s$", re_privateuse)
        ind <- grepl(pat, x, perl = TRUE)
        if(any(ind)) {
            y[pos[ind]] <-
                as.list(sprintf("Privateuse=%s",
                                substring(x[ind], 3L)))
            x[ind] <- ""
            pos <- pos[!ind]
        }
    }

    ## Now for the real thing.

    ## Remaining tags should now be as follows:
    ##   (language
    ##    ["-" script]
    ##    ["-" region]
    ##    *(["-" variant])
    ##    *(["-" extension])
    ##    ["-" privateuse]
    ## where
    ##   language    = (2*3ALPHA [-extlang])  ; shortest ISO 639 code
    ##                  / 4ALPHA              ; reserved for future use
    ##                  / 5*8ALPHA            ; registered language subtag
    ##   extlang     = *3("-" 3*ALPHA)        ; reserved for future use
    ##   script      = 4ALPHA                 ; ISO 15924 code
    ##   region      = 2ALPHA                 ; ISO 3166 code
    ##                 / 3DIGIT               ; UN M.49 code
    ##   variant     = 5*8alphanum            ; registered variants
    ##                 / (DIGIT 3alphanum)
    ##   extension   = singleton 1*("-" (2*8alphanum))
    ##   singleton   = %x41-57 / %x59-5A / %x61-77 / %x79-7A / DIGIT
    ##               ; "a"-"w" / "y"-"z" / "A"-"W" / "Y"-"Z" / "0"-"9"

    ## We handle language/extlang a bit differently (more generously).

    re_extlang <- "[[:alpha:]]{3}"
    re_language <-
        sprintf("[[:alpha:]]{2,3}(-%s){0,3}|[[:alpha:]]{4,8}", re_extlang)
    re_script <- "[[:alpha:]]{4}"
    re_region <- "[[:alpha:]]{2}|[[:digit:]]{3}"
    re_variant <- "[[:alnum:]]{5,8}|[[:digit:]][[:alnum:]]{3}"
    re_singleton <- "[abcdefghijklmnopqrstuvwyzABCDEFGHIJKLMNOPQRSTUVWYZ0123456789]"
    re_extension <- sprintf("(%s)(-[[:alnum:]]{2,8}){1,}", re_singleton)

    bad <- integer()

    if(length(pos)) {
        pat <- sprintf("^(%s)(-.*|$)", re_language)
        ind <- grepl(pat, x, perl = TRUE)
        if(!all(ind)) {
            bad <- which(!ind)
            x[bad] <- ""
        }
        y[pos[ind]] <-
            lapply(strsplit(sub(pat, "\\1", x[ind], perl = TRUE),
                            "-", fixed = TRUE),
                   function(e) {
                       c(sprintf("Language=%s", e[1L]),
                         sprintf("Extension=%s", e[-1L]))
                   })
        x[ind] <- sub(pat, "\\3", x[ind], perl = TRUE)
        ind <- nzchar(x)
        pos <- pos[ind]
        x <- x[ind]
    }

    if(length(pos)) {
        repeat {
            ## Use a loop so that we can finally stop when done.

            ## Script.
            pat <- sprintf("^-(%s)(-.*|$)", re_script)
            if(any(ind <- grepl(pat, x, perl = TRUE))) {
                y[pos[ind]] <-
                    Map(c,
                        y[pos[ind]],
                        sprintf("Script=%s",
                                sub(pat, "\\1", x[ind], perl = TRUE)))
                x[ind] <- sub(pat, "\\2", x[ind], perl = TRUE)
                ind <- nzchar(x)
                pos <- pos[ind]
                x <- x[ind]
                if(!length(x)) break
            }
        
            ## Region.
            pat <- sprintf("^-(%s)(-.*|$)", re_region)
            if(any(ind <- grepl(pat, x, perl = TRUE))) {
                y[pos[ind]] <-
                    Map(c,
                        y[pos[ind]],
                        sprintf("Region=%s",
                                sub(pat, "\\1", x[ind], perl = TRUE)))
                x[ind] <- sub(pat, "\\2", x[ind], perl = TRUE)
                ind <- nzchar(x)
                pos <- pos[ind]
                x <- x[ind]
                if(!length(x)) break
            }
        
            ## Variant(s).
            pat <- sprintf("^-(%s)(-.*|$)", re_variant)
            while(any(ind <- grepl(pat, x, perl = TRUE))) {
                y[pos[ind]] <-
                    Map(c,
                        y[pos[ind]],
                        sprintf("Variant=%s",
                                sub(pat, "\\1", x[ind], perl = TRUE)))
                x[ind] <- sub(pat, "\\2", x[ind], perl = TRUE)
                ind <- nzchar(x)
                pos <- pos[ind]
                x <- x[ind]
            }
            if(!length(x)) break

            ## Extension(s).
            pat <- sprintf("^-%s(-.*|$)", re_extension)
            while(any(ind <- grepl(pat, x, perl = TRUE))) {
                ## <NOTE>
                ## We keep the singleton prefix: this could be used in
                ## expansions of registered extensions: currently,
                ##   BCP 47 Extension U <http://tools.ietf.org/html/rfc6067>
                ##   BCP 47 Extension T <http://tools.ietf.org/html/rfc6497>
                y[pos[ind]] <-
                    Map(c,
                        y[pos[ind]],
                        sprintf("Extension=%s",
                                sub(pat, "\\1\\2", x[ind], perl = TRUE)))
                ## </NOTE>
                x[ind] <- sub(pat, "\\3", x[ind], perl = TRUE)
                ind <- nzchar(x)
                pos <- pos[ind]
                x <- x[ind]
            }
            if(!length(x)) break

            ## Private use.
            pat <- sprintf("^-%s(-.*|$)", re_privateuse)
            if(any(ind <- grepl(pat, x, perl = TRUE))) {
                y[pos[ind]] <-
                    Map(c,
                        y[pos[ind]],
                        sprintf("Privateuse=%s",
                                substring(sub(pat, "\\1", x[ind],
                                              perl = TRUE),
                                          2L)))
                x[ind] <- sub(pat, "\\4", x[ind], perl = TRUE)
            }

            break
        }
    }

    ## Be a nuisance: singletons for extensions must not be duplicated.
    ind <- as.logical(lapply(y, function(e) {
        e <- grep("^Extension=", e, value = TRUE)
        if(!length(e)) return(FALSE)
        any(duplicated(sub("^Extension=(.).*", "\\1", e)))
    }))
    if(any(ind))
        bad <- c(bad, which(ind))
    if(any(ind <- nzchar(x))) {
        bad <- c(bad, pos[ind])
    }
    if(length(bad)) {
        stop("Invalid language tag(s):",
             paste("\n ", names(y)[bad], collapse = " "),
             call. = FALSE)
    }

    if(!expand) return(y)

    x <- tolower(unlist(y))
    pos <- match(x, IANA_language_subtag_registry$Index)
    z <- IANA_language_subtag_registry$Description[pos]
    ## Special case private use ranges.
    if(!all(sapply(z, length))) {
        pos <- match(x,
                     IANA_language_subtag_registry_private_use_index_table)
        z[pos > 0L | grepl("^privateuse=", x)] <- "Private use"
    }
    z <- Map(`names<-`,
             split(z, rep.int(seq_along(y), sapply(y, length))),
             y)
    names(z) <- names(y)
    z
}

get_IANA_language_subtag_registry <-
function(con = "http://www.iana.org/assignments/language-subtag-registry")
{
    ## This is a collection of records in tag-value format, but
    ## unfortunately separated by '%%' lines rather than empty lines, so
    ## we cannot use read.dcf() directly.  Let us keep things simple:
    ## extract the records, write them out as DCF, and call read.dcf().

    lines <- readLines(con)
    ## The first line is something like
    ##   File-Date: 2009-03-13
    ## which we drop for reading the records.
    fdate <- sub(".*: *", "", lines[1L])
    pos <- grep("^%%", lines)
    lines[c(seq_len(pos[1L]), pos[-1L])] <- ""
    tcon <- textConnection(lines, encoding = "UTF-8")
    on.exit(close(tcon))
    db <- read.dcf(tcon, all = TRUE)
    ## Add index for lookups.
    subtag <- db$Subtag
    db$Index <-
        tolower(sprintf("%s=%s",
                        db$Type,
                        ifelse(is.na(subtag), db$Tag, subtag)))
    db$Type <- factor(db$Type)
    attr(db, "File_Date") <- fdate
    db
}

IANA_language_subtag_registry_language_private_use_subtags <-
    outer(letters[1L : 20L], letters,
          function(u, v) sprintf("q%s%s", u, v))
IANA_language_subtag_registry_script_private_use_subtags <-
    outer(c("a", "b"), letters[1L : 24L],
          function(u, v) sprintf("Qa%s%s", u, v))
IANA_language_subtag_registry_region_private_use_subtags <-
    c(sprintf("Q%s", LETTERS[13L : 26L]),
      sprintf("X%s", LETTERS))

IANA_language_subtag_registry_private_use_index_table <-
    tolower(c(sprintf("Language=%s",
                      IANA_language_subtag_registry_language_private_use_subtags),
              sprintf("Script=%s",
                      IANA_language_subtag_registry_script_private_use_subtags),
              sprintf("Region=%s",
                      IANA_language_subtag_registry_region_private_use_subtags)))

