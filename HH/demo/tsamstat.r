## Examples using Time Series Plotting Functions

## Displays for Direct Comparison of ARIMA Models
##
## The American Statistician, May 2002, Vol. 56, No. 2, pp. 131-138
## Richard M. Heiberger, Temple University
## Paulo Teles, Faculdade de Economia do Porto

## Many comments, marked HT, in this file refer to figure numbers in
## the Heiberger and Teles article.  Main titles for figures in the
## article have been suppressed for HH, but are left here as comments.

## For work on the computer screen, we use the residplot (called by
## tsdiagplot) and seqplot (called by tsacfplots) default lwd=0 and
## then export the graph to PostScript with the export.eps() command.
## For good quality reproduction for the book, we needed to use the
## optional lwd=1 and write directly to the PostScript file.

## co2 in S-Plus is monthly data from 1959 through 1990
## co2 in R is monthly data from 1959 through 1997
## The two co2 datasets are similar over the overlapping time periods

if.R(r=data(co2), ## Use the builtin R data set
     s={
       co2 <- as.rts(co2)                      # as.rts() is necessary
       attr(attr(co2, "tspar"), "units") <- "months"
     })

##HT New 1 1a 2
##HT	   main="Figure 1.  co2"
## export.eps(hh("tser/figure/tsamsta1.ps"))
tsacfplots(co2,                   xlab=NULL, ylab=list("co2", cex=1.5),                   lag.max=36, cex=.75, main="", lwd=.5        )
## export.eps(hh("tser/figure/tsamsta1.color.ps"))

##HT main="Figure 1a (not in printed paper).  diff(co2,1)"
## export.eps(hh("tser/figure/tsamsta1a.ps"))
tsacfplots(diff(co2,1),           xlab=NULL, ylab=list("diff(co2,1)", cex=1.5),           lag.max=36, cex=.75, main="", lwd=.5        )
## export.eps(hh("tser/figure/tsamsta1a.color.ps"))

##HT main="Figure 2.  diff(diff(co2,1), 12)"
## export.eps(hh("tser/figure/tsamsta2.ps"))
tsacfplots(diff(diff(co2,1), 12), xlab=NULL, ylab=list("diff(diff(co2,1), 12)", cex=1.5), lag.max=36, cex=.75, main="", lwd=.5        )
## export.eps(hh("tser/figure/tsamsta2.color.ps"))


###HT lag.0=FALSE
##HT main="Figure 2alt (not in printed paper).  diff(diff(co2,1), 12)  with lag.0=FALSE"
tsacfplots(diff(diff(co2,1), 12), xlab=NULL, ylab=list("diff(diff(co2,1), 12)", cex=1.5), lag.max=36, cex=.75, main="", lwd=.5, lag.0=FALSE)


##HT table 0 not in paper
cat("Table 0\n")
ddco2.loop <- if.R(s=
                   arma.loop(co2,
                             list(list(order=c(2,1,2)),
                                  list(order=c(0,1,0), period=12)))
                   ,r=
                   arma.loop(co2,
                             order=c(2,1,2),
                             seasonal=list(order=c(0,1,0), period=12))
                   )
print(ddco2.loop)


##HT 3
##HT main="Figure 3"
ddco2.diags <- diag.arma.loop(ddco2.loop, co2, lag.max=36)
ddco2.diagplot <-
tsdiagplot(armas=ddco2.loop, diags=ddco2.diags,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           main="", lwd=1)
print(ddco2.diagplot)   ## strip labels in R need 14in wide window.
## export.eps(hh("tser/figure/tsamsta3.ps"))


##HT simplest specification of Figure 4.
## this has lots of overprinting of the axis labels.
if.R(s=
     tsdiagplot(co2, model=
                list(list(order=c(2,1,2)),
                     list(order=c(0,1,1), period=12)),
                main="Figure 4simple (not in printed paper)  with defaults")
     ,r=
     tsdiagplot(co2,   ## strip labels in R need 14in wide window.
                order=c(2,1,2),
                seasonal=list(order=c(0,1,1), period=12),
                main="Figure 4simple (not in printed paper)  with defaults")
     )

## table 1
ddco2.loopPQ <-
  if.R(s=
       arma.loop(co2, list(list(order=c(2,1,2)),
                           list(order=c(0,1,1), period=12)))
       ,r=
       arma.loop(co2,
                 order=c(2,1,2),
                 seasonal=list(order=c(0,1,1), period=12))
       )
cat("Table 1\n")
print(ddco2.loopPQ)

ddco2.diagsPQ <- diag.arma.loop(ddco2.loopPQ, co2, lag.max=36)
##HT main="Figure 4"
ddco2.diagplot <-
tsdiagplot(armas=ddco2.loopPQ, diags=ddco2.diagsPQ,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           main="", lwd=1)
ddco2.diagplot  ## strip labels in R need 14in wide window.
## export.eps(hh("tser/figure/tsamsta4.ps"))

ddco2.diagplot.0 <-
tsdiagplot(armas=ddco2.loopPQ, diags=ddco2.diagsPQ,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           lag.0=FALSE,
           main="Figure 4alt (not in printed paper)   with lag.0=FALSE.")
print(ddco2.diagplot.0)  ## strip labels in R need 14in wide window.


## table 1a
cat("Table 1a (not in printed paper).\n")
print(ddco2.loopPQ[["1","1"]])

## table 2
cat("Table 2\n")
co2.arima <- ddco2.loopPQ[["0","1"]]
##The final model is equivalent to
##co2.arima <-
##  arima.mle(co2, model=list(list(order=c(0,1,1)),
##                            list(order=c(0,1,1), period=12)))
print(co2.arima)

##forecast
co2.forecast <-
  if.R(s=
       arima.forecast(co2, model=co2.arima$model, n=12)
       ,r=
       predict(co2.arima, n.ahead = 12, se.fit = TRUE)
       )

##5
##HT main="Figure 5  co2 --- 1990 observed, 1991 forecast + 95% CI"
co2.plot.forecast <-
  if.R(r={
    co2.last.year <- ts(co2[457:468],
                        start=time(co2)[457],
                        frequency=frequency(co2))
    seqplotForecast(co2.last.year, co2.forecast,
                    x.at=seq(1997,1999,.5),
                    x.labels=c(1997,"",1998,"",1999),
                    xlim=c(1996.9,1999),
                    ylim=c(360, 370),
                    main="", lwd=1, cex=2)
  }, s={
    seqplotForecast(co2[373:384], co2.forecast,
                    x.at=seq(1990,1992,.5),
                    x.labels=c(1990,"",1991,"",1992),
                    xlim=c(1990,1992),
                    main="", lwd=1, cex=1.5)
  })
## trellis.device(file=hh("tser/figure/tsamsta5.ps"), postscript, horizontal=TRUE); strip.background0()
co2.plot.forecast
## dev.off()
## export.eps(hh("tser/figure/tsamsta5.ps"))
