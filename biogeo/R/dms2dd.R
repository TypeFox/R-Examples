dms2dd <-
function (dd, mm, ss, ns) 
{
    n <- match(toupper(ns), LETTERS)
    we1 <- n == 19 | n == 23
    sgn <- ifelse(we1 == T, -1, 1)
    mm1 <- ifelse(is.na(mm), 0, mm)
    ss1 <- ifelse(is.na(ss), 0, ss)
    decdeg <- (dd + ((mm1 * 60) + (ss1))/3600) * sgn
}
