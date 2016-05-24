###
### $Id: fileparts.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Return filename parts.
###


##-----------------------------------------------------------------------------
fileparts <- function(pathname) {
    if (!is.character(pathname)) {
        stop(sprintf("argument %s must be character", sQuote("pathname")))
    } else if (!(length(pathname) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("pathname")))
    }

    ## R only expands a single tilde, optionally followed by sep
    wouldExpandTilde <- function(pathname) {
        possible <- substr(pathname, 1, 2)
        ((nchar(possible) == 1 && possible == "~") ||
         (nchar(possible) == 2 && substr(possible, 2, 2) == "/"))
    }

    tildeUser <- character(0)
    hasTilde <- substr(pathname, 1, 1) == "~"
    if (hasTilde && wouldExpandTilde(pathname)) {
        ## Augment tilde with bogus value prevent expansion by path.expand
        luser <- "xxxxxx"    # :HACK: assumed not to exist
        tildeUser <- paste("~", luser, sep = "")
        pathname <- sub("~", tildeUser, pathname)
    }

    fname <- basename(pathname)
    if (fname == pathname) {
        pathstr <- ""
    } else {
        hasTrailingSep <- function(pathname) {
            nchars <- nchar(pathname)
            lastChar <- substr(pathname, nchars, nchars)
            ## :TODO: There must be a standard function for the following...
            if (.Platform$OS.type == "windows") {
                lastChar == "/" || lastChar == "\\"
            } else {
                lastChar == "/"
            }
        }

        pathstr <- if (hasTrailingSep(pathname)) {
                       fname <- ""
                       ## dirname normalizes path ending with fsep so
                       ## append a character so trailing sep is kept
                       dirname(paste(pathname, "x", sep = ""))
                   } else {
                       dirname(pathname)
                   }
    }

    if (fname == ".") {
        ## Handle relative current directory
        name <- character(0)
        ext <- "."
    } else if (fname == "..") {
        ## Handle relative parent directory
        name <- "."
        ext <- "."
    } else {
        name <- {
                    ext.re <- "\\.[^\\.]*$" # match from last period to end
                    sub(ext.re, "", fname)
                }
        if (name == "") {
            ## Handle UNIX hidden files
            name <- character(0)
            ext <- fname
        } else {
            split.re <- name
            ext <- unlist(strsplit(fname, split.re, fixed = TRUE))[2]
        }
    }

    if (length(tildeUser) > 0) {
        switch(EXPR = charmatch("~", c(pathstr, name)),
               pathstr <- sub(tildeUser, "~", pathstr),
               name <- sub(tildeUser, "~", name))
    }

    return(list(pathstr = pathstr,
                name    = name,
                ext     = ifelse(!is.na(ext), ext, ""),
                versn   = ""))
}

