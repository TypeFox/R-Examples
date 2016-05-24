dms2rad <-
function(d = '33d 09m 35.0s') {
# Conversion from degrees, minutes second string to radians
# Input must be string for case when degree term is between 0 and 1.
# Sign character may be W or E for longitude, also P or M for plus/minus.
#
# A. Harris

    # Strip out blanks from string, substitute for W and E if longitude
    d <- gsub(' ','', d)
    d <- gsub('[Ww]','-', d)
    d <- gsub('[Ee]','', d)
    # Get length of input vector
    len <- length(d)
    rad <- numeric(len)
    # Loop over elements
    for (i in 1:len) {
        sgn <- ifelse(substr(d[i], 1, 1)=='-', -1, 1)
        # Convert to numeric, dividing on any of h, d, m, s, colon, or comma
        dms <- as.numeric(strsplit(d[i], '[hdms:,]')[[1]])
        # Fill out array with zeros as needed
        dms <- c(dms, numeric(3-length(dms)))
        # Compute absolute radians
        rad[i] <- (abs(dms[1]) + (abs(dms[2]) + abs(dms[3])/60)/60)*pi/180
        # Return signed radians
        rad[i] <- sgn*rad[i]
    }
    return(rad)
}
