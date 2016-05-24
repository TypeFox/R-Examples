syms <-
function (z, zrange = c(NA, NA), p = 1) 
{
    dmin <- 0.1
    dmax <- 1
    ddiff <- dmax - dmin
    if (is.na(zrange[1])) 
        zmin <- min(z)
    else zmin <- zrange[1]
    if (is.na(zrange[2])) 
        zmax <- max(z)
    else zmax <- zrange[2]
    zdiff <- zmax - zmin
    z <- ifelse(z < zmin, zmin, z)
    z <- ifelse(z > zmax, zmax, z)
    zdiam <- dmin + ddiff * ((z - zmin)/zdiff)^p
    invisible(zdiam)
}
