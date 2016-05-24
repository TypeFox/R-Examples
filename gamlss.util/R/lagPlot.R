# created 6-1-12
# to do 
# when plotting use the actual variable names
lagPlot <- function(y, x=NULL, lags=0, corr = TRUE, smooth = TRUE)
{
lag.plot1 <- function (series, lags = 1, corr = TRUE, smooth = TRUE) 
{
    name1 = paste(deparse(substitute(series)), "(t-", sep = "")
    name2 = paste(deparse(substitute(series)), "(t)", sep = "")
    data1 = as.ts(series)
    lags = as.integer(lags)
    prow = ceiling(sqrt(lags))
    pcol = ceiling(lags/prow)
    a = acf(series, lags, plot = FALSE)$acf[-1]
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(prow, pcol), mar = c(2.5, 4, 2.5, 1), cex.main = 1.1, 
        font.main = 1)
    for (h in 1:lags) {
        plot(lag(series, -h), data1, xy.labels = FALSE, main = paste(name1, 
            h, ")", sep = ""), ylab = name2, xlab = "")
        if (smooth == TRUE) 
            lines(lowess(ts.intersect(lag(series, -h), series)[, 
                1], ts.intersect(lag(series, -h), series)[, 2]), 
                col = "red")
        if (corr == TRUE) 
            legend("topright", legend = round(a[h], digits = 2), 
                text.col = "blue", bg = "white", x.intersp = 0)
        on.exit(par(old.par))
    }
}

lag.plot2<-function (series1, series2, lags = 0, corr = TRUE, smooth = TRUE) 
{
    name1 = paste(deparse(substitute(series1)), "(t-", sep = "")
    name2 = paste(deparse(substitute(series2)), "(t)", sep = "")
    data1 = as.ts(series1)
    data2 = as.ts(series2)
    lags = as.integer(lags)
    m1 = lags + 1
    prow = ceiling(sqrt(m1))
    pcol = ceiling(m1/prow)
    a = ccf(series1, series2, lags, plot = FALSE)$acf
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(prow, pcol), mar = c(2.5, 4, 2.5, 1), cex.main = 1.1, 
        font.main = 1)
    for (h in 0:lags) {
        plot(lag(series1, -h), series2, xy.labels = FALSE, main = paste(name1, 
            h, ")", sep = ""), ylab = name2, xlab = "")
        if (smooth == TRUE) 
            lines(lowess(ts.intersect(lag(series1, -h), series2)[, 
                1], ts.intersect(lag(series1, -h), series2)[, 
                2]), col = "red")
        if (corr == TRUE) 
            legend("topright", legend = round(a[m1 - h], digits = 2), 
                text.col = "blue", bg = "white", x.intersp = 0)
        on.exit(par(old.par))
    }
}
  
  
if (is.null(x))
   {lags <- if (lags==0) 1  else lags
   lag.plot1(series=y, lags=lags,  corr = corr, smooth = corr)
   }
  else
   lag.plot2(series1=x, series2=y, lags=lags,  corr = corr, smooth = corr)  
}
