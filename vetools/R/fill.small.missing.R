# Verified 1.3.18
## fill.small.missing (vetools)
## 18 Nov. 2013
## A. M. Sajo-Castelli
## D. Villalta
## Version 2.0
## Minor revisions and window bug fixed.
fill.small.missing<-function (serie, max.len = 3 * 30, func = median) 
{
        if ( frequency(serie) != 365.25 ) stop ("Series must be of frequency 365.25.")
        if (any(serie[!is.na(serie)] < 0)) {
                stop("Series must be non negative.")
        }
        time.i = time(serie)
        table.ref = (serie >= 0)
        table.ref[is.na(table.ref)] = FALSE
        END = length(table.ref)
        idx_s = 1
        AA = seq(floor(start(serie)[1]), floor(end(serie)[1]), 1)
        while (idx_s < END) {
                idx_s = idx_s - 1 + match(FALSE, table.ref[idx_s:END])
                if (is.na(idx_s)) {
                        idx_e = END
                        break
                }
                idx_e = (idx_s) + match(TRUE, table.ref[(idx_s + 1):END]) - 1
                if (is.na(idx_e)) {
                        idx_e = idx_s
                        break
                }
                if ((idx_e - idx_s + 1) > max.len) {
                        warning(paste("Missing gap grater than", max.len, 
                                      "days, skipping..."))
                        idx_s = idx_e + 1
                        next
                }
                try(length(window(serie, start = time.i[idx_s], end = time.i[idx_e])) != (idx_e - idx_s + 1))
                Mediana = array(numeric(0), dim = c(length(AA), (idx_e - idx_s + 1)))
                conty = 0
                for (y in AA) {
                        conty = conty + 1
                        if (y == tis::year(time.i[idx_s])) {
                                Mediana[conty, ] = NA
                                next
                        }
                        lstart = time.i[idx_s] - floor(time.i[idx_s]) + y - 1 /  (2 * 365.25)
                        lend = time.i[idx_e] - floor(time.i[idx_s]) + y + 1 / (2 * 365.25)
                        if ( lstart > time.i[length(time.i)] ) { break }
                        w = window(serie, start = lstart, end = lend, extend = T)
                        Mediana[conty, ] = w[1:(idx_e - idx_s + 1)]
                }
                Med = apply(Mediana, 2, func, na.rm = T)
                serie[idx_s:idx_e] <- Med
                idx_s = idx_e + 1
        }
        serie = window(serie, start = time.i[1], end = time.i[idx_e])
        return(serie)
}