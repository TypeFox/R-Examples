dmsparsefmt <-
function (x, fmt) 
{
    z <- {
    }
    for (i in 1:length(x)) {
        z1 <- dmsabs(x[i], fmt)
        mm <- ifelse(is.na(z1$mm), 0, z1$mm)
        ss <- ifelse(is.na(z1$ss), 0, z1$ss)
        decdeg <- dms2dd(z1$dd, mm, ss, z1$L)
        z2 <- data.frame(z1, decdeg, stringsAsFactors = F)
        z <- rbind(z, z2)
    }
    z <- data.frame(coord = x, z, stringsAsFactors = F)
    return(z)
}
