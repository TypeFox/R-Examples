date2days <- function(date, dateformat, startmonth)
{
    days <- numeric()
    ndays <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    posmm <- regexpr("mm", as.character(dateformat), ignore.case = "TRUE")[1]
    posdd <- regexpr("dd", as.character(dateformat), ignore.case = "TRUE")[1]
    for (n in 1:length(date)) {
        mm <- as.numeric(substr(date[n], posmm, posmm + 1))
        dd <- as.numeric(substr(date[n], posdd, posdd + 1))
        d <- dd
        if (startmonth < mm) {
            for (i in startmonth:(mm - 1)) {
                d <- d + ndays[i]
            }
        }
        if (startmonth > mm) {
            for (i in (startmonth):12) {
                d <- d + ndays[i]
            }
            if (mm > 1) {
                for (i in 1:(mm - 1)) {
                  d <- d + ndays[i]
                }
            }
        }
        days <- c(days, d)
    }
    return(days)
}