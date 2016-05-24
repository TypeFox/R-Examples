uniqueformats <-
function (x2) 
{
    fm <- getformat(x2)
    ufm <- unique(fm)
    dat <- {
    }
    for (i in 1:length(ufm)) {
        f <- which(fm == ufm[i])
        f1 <- f[1]
        d <- data.frame(coord = x2[f1], format = fm[f1], stringsAsFactors = F)
        dat <- rbind(dat, d)
    }
    return(dat)
}
