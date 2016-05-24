j2b <-
function(ra='17:30:30', dec='-28:47'){
# Convert from J2000 to B1950 coordinates
# Input RA and Dec as strings, elements separated by commas
#
# From the Explanatory supplement to the Astronomical
# Almanac, Seidelmann (ed.), c 1992, chapter 3.213
#
# A. Harris 2010.2.11, 2012.6.24

    # Requires hms2rad, rad2hms, dms2rad, rad2dms

    # Dates: t=-0.5 (50 years before epoch 2000.0), T=0 (epoch 2000.0) (3.211-2)
    t = -0.5

    # Scaling factors M and N (3.213-2), in radians:
    M <- (1.2812323 + (0.0003879 + 0.0000101*t)*t)*t*pi/180
    N <- (0.5567530 - (0.0001185 - 0.0000116*t)*t)*t*pi/180

    # Reduction from J2000.0 (0 subscript) (2.213-1)
    alpha.0 <- hms2rad(ra)
    delta.0 <- dms2rad(dec)
    alpha.m <- alpha.0 + (M + N*sin(alpha.0)*tan(delta.0))/2
    delta.m <- delta.0 + N*cos(alpha.m)/2
    alpha <- alpha.0 + M + N*sin(alpha.m)*tan(delta.m)
    delta <- delta.0 + N*cos(alpha.m)

    # Return B1950 RA and Dec, round seconds for reasonable precision
    out <- list(ra1950=rad2hms(alpha, 1), dec1950=rad2dms(delta, 0))
    class(out) <- 'j2b'
    out
}

# Print method for j2b class
# AH 2012.05.20

print.j2b <- function(x, digits = NULL, quote = TRUE, na.print = NULL,
                          print.gap = NULL, right = FALSE, max = NULL,
                          useSource = TRUE, ...) {
    cat(x$ra1950, x$dec1950, '\n')
}
