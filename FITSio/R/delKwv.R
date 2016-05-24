delKwv <-
function(keyw, headerName) {
# Function deletes 'keyword = value /note' card image from FITS header
#
# For multi-card select, use grep, e.g. idx <- grep('^TEST[7-9] ', header)
#
# A. Harris 2012.10.13

    # find keyword; '=' in col. 9 defines as keyword
    keyw <- sprintf('%-8s=', substr(keyw, 1, 8))
    idx <- which(toupper(keyw) == substr(toupper(headerName), 1, 9))
    if (length(idx) == 0) stop('*** No such keyword to delete ***')

    # eliminate card image
    headerName <- headerName[-idx]

    # return modified structure
    headerName
}
