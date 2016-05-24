parseHdr <-
function(headerName) {
## Take header card images and produce a vector with elements keyword, value
##
## A. Harris 2012.10.30, 2013.01.28

    # identify keyword=value pairs by checking for '= ' as 9th and 10th
    # characters (strict FITS encoding)
    idx <- which(substr(headerName, 9, 10)=='= ')

    # elminate notes (/ is not a valid keyword character)
    for (i in idx) {
        headerName[i] <- strsplit(headerName[i], '/')[[1]][1]
    }
    idx <- 1:length(headerName)

    # parse headers
    hdr <- unlist(strsplit(headerName[idx], '='))

    # extract string keywords
    idx <- grep("'", hdr)
    for (i in idx) {
        # replace doubled single quotes with improbable dummy string
        hdr[i] <- gsub("''", 'aAlJ2fZ47xx', hdr[i])
        # extract keyword
        hdr[i] <- strsplit(hdr[i], "'")[[1]][2]
        # replace improbable dummy string with doubled single quotes
        hdr[i] <- gsub('aAlJ2fZ47xx', "''", hdr[i])
    }

    # eliminate leading and trailing spaces
    for (i in 1:length(hdr)) {
        hdr[i] <- sub('^ *', '', hdr[i])
        hdr[i] <- sub(' *$', '', hdr[i])
    }

    # return parsed header
    hdr
}
