addKwv <-
function(keyw, val, note='', headerName) {
# Function adds 'keyword = value /note' card image to FITS header
#
# A. Harris 2012.10.13
    # check for valid inputs, stop otherwise
    if (!is.character(keyw)) stop('*** Keyword must be a string ***')
    if (!is.character(note)) stop('*** Note must be a string ***')
    # format keyword and note
    keyw <- toupper(strtrim(keyw, 8))
    note <- strtrim(note, 47)

    # modify card images, header, and header index
    headerName<- c(headerName, newKwv(keyw, val, note))

    # return modified structure
    headerName
}
