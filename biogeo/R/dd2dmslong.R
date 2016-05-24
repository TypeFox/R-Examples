dd2dmslong <-
function (decdeg) 
{
    dd <- floor(abs(decdeg))
    m1 <- (abs(decdeg) - dd) * 60
    mm <- floor(m1)
    ss <- round((m1 - floor(m1)) * 60, 1)
    for (i in 1:length(ss)) {
        if (ss[i] >= 60) {
            mm[i] <- mm[i] + 1
            ss[i] <- 0
        }
    }
    for (i in 1:length(mm)) {
        if (mm[i] >= 60) {
            dd[i] <- dd[i] + 1
            mm[i] <- 0
        }
    }
    ns <- ifelse(abs(decdeg) - decdeg == 0, "E", "W")
    dms <- data.frame(dd, mm, ss, ns)
    names(dms) <- c("xdeg", "xmin", "xsec", "EW")
    return(dms)
}
