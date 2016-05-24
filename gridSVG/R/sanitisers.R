#
# Takes a grob name as an input and returns a selector capable of being
# used as a CSS selector in JS, e.g.:
# jQuery - $(...)
# D3 - d3.select(...)
# document.querySelector(...)
#
# Note though, document.getElementById() works *without* escaping.
#
# Mathias Bynens has written a *JS* implementation of this escaping which we
# could use for any JS code, but I think an R port is ideal. This is because
# we are not dependent on any external library, and the data is fixed. This
# means that it can easily be abstracted away for most purposes.
#
# https://github.com/mathiasbynens/mothereff.in/blob/master/css-escapes/eff.js
#
# No escaping necessary for XPath! Should just be "//svg:*[@id='OUR_ID']"
# Only need to escape quoting character, (') by default
#
# Params:
# x, a grob name to escape
# escapeNonASCII, if we have unicode characters, should we escape them
# escapeJS, if we want to have safe values for JS, we need to escape
#           that too. This is FALSE by default because RJSONIO does this
#           for us.
escapeSelector <- function(x, escapeNonASCII = FALSE, escapeJS = FALSE) {
    isUTF8 <- localeToCharset()[1] == "UTF-8"
    # *must* do \ first because it is the escape character and we do not
    # want to escape any escape characters that we would have added in.

    escapeCache <- list(
        "\b" = "\\b",
        "\t" = "\\t",
        "\n" = "\\n",
        "\v" = "\\x0b",
        "\f" = "\\f",
        "\r" = "\\r",
        "\\" = "\\\\",
        "'" = "\\'",
        '"' = '\\"'
    )

    # If we have a unicode string, we have to do some special checking
    if (isUTF8) {
        # Line separator
        escapeCache[["\u2028"]] <- "\\u2028"
        # Paragraph separator
        escapeCache[["\u2029"]] <- "\\u2029"
    }

    chars <- substring(x, 1:nchar(x), 1:nchar(x))
    n <- length(chars)
    # Begin CSS escaping
    for (i in 1:n) {
        charcode <- as.integer(charToRaw(chars[i]))
        # Some unicode char or non-printable char, just escape it
        if (escapeNonASCII && (length(charcode) > 1 || charcode < 32 || charcode > 126))
            chars[i] <- paste0("\\", chars[i])
        else {
            if (grepl("[\t\n\v\f:]", chars[i])) {
                chars[i] <- paste0("\\", paste0(as.hexmode(charcode), collapse = ""), " ")
            } else if (grepl("[ !\"#$%&'()*+,./;<=>?@\\[\\\\\\]^`\\{|\\}~]", chars[i], perl = TRUE)) {
                chars[i] <- paste0("\\", chars[i])
            } else {
                # Do nothing, no escaping required
            }
        }
    }

    chars <- paste0(chars, collapse = "")

    # We shouldn't need to escape anything because RJSONIO does this
    # by default, but the code is here if required
    if (! escapeJS)
        return(paste0("#", chars))

    chars <- substring(chars, 1:nchar(chars), 1:nchar(chars))
    n <- length(chars)
    # Begin JS escaping
    for (i in 1:n) {
        # Would ideally do grepl([\x20-\x26\x28-\x5b\x5d-\x7e], ...)
        # In decimal this is: 32-38 40-91 93-126
        matchVector <- c(32:38, 40:91, 93:126)
        charcode <- as.integer(charToRaw(chars[i]))
        if (length(charcode) == 1 && charcode %in% matchVector) {
            # Do nothing, we do not need to escape this character
        } else if (! is.null(escapeCache[[chars[i]]])) {
            chars[i] <- escapeCache[[chars[i]]]
        } else {
            chars[i] <- escapeCache[[chars[i]]] <- paste0("\\", chars[i])
        }
    }

    paste0("#", paste0(chars, collapse = ""))
}

# Much simpler for XPath, just escape quote chars
# x, a grob name to escape
escapeXPath <- function(x) {
    # We're going to use single quotes, so escape only those.
    # RJSONIO uses ", so is safer to use single quotes
    if (grepl("'", x))
        x <- gsub("'", "\\\\'", x)

    # Look through the entire document for all elements,
    # then return only the element with the matching ID.
    # This could possibly collide with existing IDs if
    # embedded within a document (e.g. HTML) but very unlikely.
    paste0("//*[@id='", x, "']")
}

