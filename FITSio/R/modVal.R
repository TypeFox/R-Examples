modVal <-
function(keyw, val, note='', headerName) {
# Function modifies value in 'keyword = value /note' card image from FITS header
#    and allows overwrite of existing note
#
# A. Harris 2012.10.13

    # keyword for all card images
    keyw <- sprintf('%-8s=', substr(keyw, 1, 8))
    idx <- which(toupper(keyw) == substr(toupper(headerName), 1, 9))
    if (length(idx) == 0) stop('*** No such keyword to modify ***')
    if (length(idx) > 1) stop('*** Multiple keywords to modify ***')

    # parse card structure
    tmp <- headerName[idx]
    tmp <- unlist(strsplit(tmp, '/'))
    if (note == '') {
        note <- ifelse(length(tmp)==2, tmp[2], '')
    } else {
        if (!is.character(note)) stop('*** Note must be a string ***')
        note <- strtrim(note, 47)
    }
    keyw <- unlist(strsplit(gsub(' ', '', tmp[1]), '='))[1]

    # write new card
    headerName[idx] <- newKwv(keyw, val, note)

    # return modified structure
    headerName
}
