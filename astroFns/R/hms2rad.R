hms2rad <-
function(h='12h 3m 45.6s') {
# Conversion from degrees, minutes second string to radians
# Input must be string for case when degree term is between 0 and 1.
#
# A. Harris

    h <- gsub(' ','', h)
    # Get length of input vector
    len <- length(h)
    rad <- numeric(len)
    # Loop over elements
    for (i in 1:len) {
        sgn <- ifelse(substr(h[i], 1, 1)=='-', -1, 1)
        # Convert to numeric, dividing on any of h, d, m, s, colon, or comma
        hms <- as.numeric(strsplit(h[i], '[hdms:,]')[[1]])
        # Fill out array with zeros as needed
        hms <- c(hms, numeric(3-length(hms)))
        # Compute absolute radians
        rad[i] <- (abs(hms[1]) + (abs(hms[2]) + abs(hms[3])/60)/60)*pi/12
        # Return signed radians
        rad[i] <- sgn*rad[i]
    }
    return(rad)
}

