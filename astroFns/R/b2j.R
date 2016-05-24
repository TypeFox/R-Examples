b2j <-
function(ra = "17h42m29.3076s", dec = "-28d59m18.484s"){
# Convert from B1950 to J2000 coordinates
#
# Input RA and Dec as strings, elements separated by commas
#
# From the Explanatory supplement to the Astronomical
# Almanac, Seidelmann (ed.), c 1992, chapter 3.213
#
# A. Harris 2010.2.11, 2012.6.12

    # Requires hms2rad, rad2hms, dms2rad, rad2dms

    # Dates: t=-0.5 (50 years before epoch 2000.0), T=0 (epoch 2000.0) (3.211-2)
    t = -0.5

    # Scaling factors M and N (3.213-2), in radians:
    M <- (1.2812323 + (0.0003879 + 0.0000101*t)*t)*t*pi/180
    N <- (0.5567530 - (0.0001185 - 0.0000116*t)*t)*t*pi/180

    # Reduction to J2000.0 (0 subscript) (2.213-1)
    alpha <- hms2rad(ra)
    delta <- dms2rad(dec)
    alpha.m <- alpha - (M + N*sin(alpha)*tan(delta))/2
    delta.m <- delta - N*cos(alpha.m)/2
    alpha.0 <- alpha - M - N*sin(alpha.m)*tan(delta.m)
    delta.0 <- delta - N*cos(alpha.m)

    # Return J2000 RA and Dec, round seconds for reasonable precision
    # Output lists
    out <- list(ra2000=rad2hms(alpha.0, 1), dec2000=rad2dms(delta.0, 0))
    class(out) <- 'b2j'
    out

}

###################################

# Print method for b2j class
# AH 2012.03.23

print.b2j <- function(x, digits = NULL, quote = TRUE, na.print = NULL,
                          print.gap = NULL, right = FALSE, max = NULL,
                          useSource = TRUE, ...) {
    cat(x$ra2000, x$dec2000, '\n')
}
