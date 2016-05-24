sep <-
function (dx) 
{
    x <- dx[, 1]
    let <- dx[, 2]
    xx <- x
    ft <- c(as.character(0:9), "[.]")
    for (i in 1:11) {
        x <- str_replace_all(x, ft[i], "0")
    }
    dat <- {
    }
    for (i in 1:length(x)) {
        fmt <- x[i]
        h <- str_locate_all(fmt, "0")
        v <- matrix(unlist(h), ncol = 2, byrow = F)
        v1 <- v[, 1]
        st <- v[1, 1]
        ed <- v[nrow(v), 1]
        del <- setdiff(st:ed, v1)
        if (length(del) == 0) {
            deg <- as.numeric(substring(xx[i], st, ed))
            min <- NA
            sec <- NA
        }
        else {
            del1 <- del[1]
            del2 <- del[length(del)]
            fd <- which(v1 < del1)
            deg <- as.numeric(substring(xx[i], fd[1], fd[length(fd)]))
            fm <- which(v1 > del1)
            fmin1 <- v1[fm[1]]
            if (del1 == del2) {
                min <- as.numeric(substring(xx[i], fmin1, ed))
                sec <- 0
                fm <- floor(min)
                if (min - fm > 0) {
                  sec <- abs(floor(min) - min) * 60
                  min <- floor(min)
                }
            }
            else {
                fm2 <- which(v1 >= fmin1 & v1 < del2)
                fm2 <- ifelse(length(fm2) > 1, fm2[length(fm2)], 
                  fm2)
                fmin2 <- v1[fm2]
                min <- as.numeric(substring(xx[i], fmin1, fmin2))
                fs <- which(v1 > fmin2)
                fsec <- v1[fs]
                sec <- as.numeric(substring(xx[i], fsec[1], ed))
            }
        }
        dt <- data.frame(deg, min, sec)
        dat <- rbind(dat, dt)
    }
    dat <- data.frame(dat, L = let, stringsAsFactors = F)
    return(dat)
}
