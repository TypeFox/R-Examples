parsecoords <-
function (dat, d1, fmtstr) 
{
    fm <- getformat(dat)
    dat1 <- data.frame(dat, fm)
    n <- length(dat)
    n1 <- nrow(d1)
    d2 <- {
    }
    ff <- {
    }
    for (i in 1:n1) {
        f <- which(dat1$fm == d1$format[i])
        s <- dat1[f, ]
        z <- dmsparsefmt(s$dat, fmtstr[i])
        d2 <- rbind(d2, z)
        ff <- c(ff, f)
    }
    d2 <- d2[ff, ]
    names(d2) <- c("dms", "deg", "min", "sec", "L", "dd")
    return(d2)
}
