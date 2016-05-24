substddmm <-
function (dc) 
{
    dd <- floor(abs(dc))
    m1 <- (abs(dc) - dd) * 60
    mm <- floor(m1)
    ss <- (m1 - floor(m1)) * 60
    ns <- ifelse(abs(dc) - dc == 0, "N", "S")
    z <- dms2dd(mm, dd, ss, ns)
}
