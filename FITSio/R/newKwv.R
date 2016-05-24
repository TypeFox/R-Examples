newKwv <-
function(keyw, val, note='') {
# Function writes new 'keyword = value /note' card image for FITS header
#
# A. Harris 2012.10.13

    # check for valid inputs, stop otherwise
    if (!is.character(keyw))
        stop('*** Non-numeric keyword must be a string ***')
    # format keyword and note
    keyw <- toupper(strtrim(keyw, 8))
    noteSep <- ifelse(note == '', ' ', '/')

    # format value; first numeric, then string (follows fixed format convention)
    if (is.numeric(val)) {
        txt <- sprintf('%-20.14g   %s %-46s', val, noteSep, note)
    } else {
        if (nchar(val) > 68) {
            val <- strtrim(val, 68)
            cat('   *** Truncated value string for keyword =', keyw,
                'to 68 characters ***\n')
            warning('Truncated value string in newKw')
        }
        txt <- sprintf('\'%-1s\' %s %-59s', val, noteSep, note)
    }
    sprintf('%-8s= %-70s', keyw, strtrim(txt, 70))
}
