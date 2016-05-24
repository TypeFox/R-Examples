readRwl = function (fname, header = NULL, n.header=NULL, info=TRUE) 
{
    dat = file(fname, "r")
    if (is.null(header)) {
        hdr1 = readLines(dat, n = 1)
        if (nchar(hdr1) < 12) 
stop("First line in .rwl file ends before col 12")
        yrcheck = suppressWarnings(as.numeric(substr(hdr1, 9, 12)))		          #FC        
        if (is.null(yrcheck) | length(yrcheck) != 1 | is.na(yrcheck) | 
            yrcheck < -10000 | yrcheck > 10000) {
            cat("There appears to be a header in the rwl file\n")
            is.head = TRUE
if (is.null(n.header)) n.header=1
        }
        else {
            is.head = FALSE
            #cat("There does not appear to be a header in the rwl file\n")
        }
        close(dat)
        dat = file(fname, "r")
    }
    else is.head = header
    if (is.head) {
        dat1 = readLines(dat, n = n.header+1)[-c(1:n.header)]                   #FC
    }
    else dat1 = readLines(dat, n = 1)
    yrcheck = as.numeric(substr(dat1, 9, 12))
    if (is.null(yrcheck) | length(yrcheck) != 1) {
        stop("Cols 9-12 of first data line not a year")
    }
    close(dat)
    skip.lines = ifelse(is.head, n.header, 0)
    dat = read.fwf(fname, c(8, 4, rep(6, 10)), skip = skip.lines, 
        strip.white = TRUE, blank.lines.skip = TRUE)
    dat = dat[!apply(is.na(dat), 1, all), ]
    series = dat[, 1]
    series.ids = unique(series)
    nseries = length(series.ids)
if (info) cat(nseries, "series\n", sep = " ")
    series.index = match(series, series.ids)
    min.year = (min(dat[, 2])%/%10) * 10
    max.year = ((max(dat[, 2]) + 10)%/%10) * 10
    span = max.year - min.year + 1
    rw.mat = matrix(NA, ncol = nseries, nrow = span)
    colnames(rw.mat) = as.character(series.ids)
    rownames(rw.mat) = min.year:max.year
    for (i in 1:nseries) {
        id = unique(dat[series.index == i, 1])
        decade.yr = dat[series.index == i, 2]
        first.yr = dat[series.index == i, 2][1]
        x = dat[series.index == i, -c(1, 2)]
        for (j in 1:nrow(x)) {
            yr = decade.yr[j]
            for (k in 1:ncol(x)) {
                if (is.na(x[j, k])) 
                  break
                rw.mat[as.character(yr), i] = x[j, k]
                yr = yr + 1
            }
        }
    }
    if (any(rw.mat == -999, na.rm = TRUE) | any(rw.mat == 999, 
        na.rm = TRUE)) {
        prec = 0.01
        rw.mat[rw.mat == 999] = NA
        rw.mat[rw.mat == -999] = NA
    }
    if (any(rw.mat == -9999, na.rm = TRUE) | any(rw.mat == 9999, 
        na.rm = TRUE)) {
        prec = 0.001
        rw.mat[rw.mat == -9999] = NA
        rw.mat[rw.mat == 9999] = NA
    }
    rw.mat = rw.mat * prec
    fix.internal.na = function(x) {
        for (i in 2:length(x)) {
            if (!is.na(x[i - 1]) & is.na(x[i]) & !is.na(x[i + 
                1])) 
                x[i] = 0
        }
        x
    }
    rw.mat = apply(rw.mat, 2, fix.internal.na)
    rw.mat = rw.mat[!apply(is.na(rw.mat), 1, all), ]
    rw.df = as.data.frame(rw.mat)
    rw.df
}
