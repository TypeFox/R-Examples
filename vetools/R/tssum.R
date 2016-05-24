# Verified 1.3.18
# Examples:
# #MENSUAL
# meses=1:12 # sea el inicio de cada temporada
# #BIMENSUAL
# meses=seq(1, 12, 2) # sea el inicio de cada temporada
# #TRIMETRAL
# meses=seq(1, 12, 3) # sea el inicio de cada temporada
# #TRIMETRAL ###
# meses=c(2, 5, 8, 11) # sea el inicio de cada temporada
# #SEMESTRAL
# meses=c(1, 7) # sea el inicio de cada temporada
# #SEMESTRAL ###
# meses=c(4, 10) # sea el inicio de cada temporada
tssum <-
function(series, months=1:12, max.na.fraction=0.3, safe.check=FALSE) {
        meses = months
        if ((attributes(series)$tsp[3] != 365.25) | (attributes(series)$tsp[3] == 365)) {stop("tssum only handles frequencies 365 and 365.25.")}
        ts2time <- function(e) {
                if (length(e) == 1) {
                        return(tis::jul(e))
                } else {
                        return(tis::jul(e[1]+e[2]/(365+lubridate::leap_year(e[1]))))
                }
        }
        per = length(meses)
        e = series
        if ( per == 1 ) {
                deltat = 12
        } else {
                deltat = meses[2] - meses[1]
        }

        M = floor(as.numeric(ceiling(difftime((ts2time(end(e))), (ts2time(start(e))), units="weeks"))/52))

        nd = rep(NA, M)
        est.men = ts(nd, start=c(tis::year(ts2time(start(e))), tis::month(ts2time(start(e)))), end=c(tis::year(ts2time(end(e))), tis::month(ts2time(end(e)))), frequency=per)
        AA = seq(tis::year(ts2time(start(e))), tis::year(ts2time(end(e))))
        loss = 0L
        for ( y in AA ) {
                co = 0
                for (p in meses){
                        p.jul = tis::jul(paste(y, p, "01", sep="-"))
                        if ( ts2time(end(series)) <= p.jul ) { break }
                        co = co + 1
                        dOY.p.jul = tis::dayOfYear(p.jul)
                        f = dOY.p.jul + diasdelmes(y, p:(p + deltat - 1)) - 1
                        dy = tis::dayOfYear(tis::jul(paste(y, "12-31", sep="-")))
                        if (f > dy ) { yf = y + 1; f = f - dy + 1 } else { yf = y }
                        w = window(e, start = c(y, dOY.p.jul), end = c(yf, f))
                        #  w = window(e, start = c(y, dOY.p.jul), end = c(yf, f), extend=T)
                        if ( (sum(is.na(w)) / length(w)) > max.na.fraction ) {
                                w.sum = NA
                                loss = loss + sum(w, na.rm=T) # Loss due to NA imputation
                        } else {
                                w.sum = sum(w, na.rm=T)
                        }
                        window(est.men, start=c(y, co), end=c(y, co)) = w.sum
                        ta1 = sum(window(e, start=start(e), c(yf, f)), na.rm=T)
                        ta2 = sum(window(est.men, start=start(est.men), end=c(y, co)), na.rm=T)
                        if ( (abs(ta1 - ta2) > 1e-3 ) & !is.na(w.sum) ) {
                                window(est.men, start=c(y, co), end=c(y, co)) = w.sum + (ta1 - ta2) - loss
                        }
                        if (safe.check) {
                                ta2 = sum(window(est.men, start=start(est.men), end=c(y, co)), na.rm=T)
                                if (abs(ta1 - ta2 - loss) > 1e-3){ stop("Mass conservation checksum error.") }
                        }
                }
        }
        return(window(est.men, start=start(series), end=end(series)))
}
