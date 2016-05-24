print.trend.test <-
function (x, ...){
##    Copyright (C) 2015, 2016  Thorsten Pohlert
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##    This is a generic print function.
##
    if (x$method =="SMK"){
        cat(paste("Seasonal Mann-Kendall Test without correlation"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("two-sided homogeinity test"),fill=TRUE)
        cat(paste("H0: S = 0 (no trend)"), fill=TRUE)
        cat(paste("HA: S != 0 (monotonic trend)"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("Statistics for individual seasons"), fill=TRUE)
        res.S <- data.frame(S = x$Sg,
                            varS = round(x$varSg,1),
                            Z = round(x$Zg,1),
                            tau = round(x$taug,3),
                            pvalue = format.pval(x$pvalg)
                            )

        print(res.S)
        cat(paste(" "), fill=TRUE)
        cat(paste("Statistics for total series"), fill=TRUE )
        res.t <- data.frame(S = x$Stot,
                            varS = x$Varianz,
                            Z = round(x$Z,1),
                            tau = round(x$tautot,3),
                            pvalue = format.pval(x$pvalue)
                            )
        print(res.t)
    }
    else if(x$method =="CSMK"){
        cat(paste("Seasonal Mann-Kendall Test with correlations"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("two-sided homogeinity test"),fill=TRUE)
        cat(paste("H0: S = 0 (no trend, seasons are correlated)"), fill=TRUE)
        cat(paste("HA: S != 0 (monotonic trend)"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
               cat(paste("Statistics for total series"), fill=TRUE )
        res.t <- data.frame(S = x$Stot,
                            varS = x$Varianz,
                            Z = round(x$Z,1),
 #                           tau = round(x$tautot,3),
                            pvalue = format.pval(x$pvalue)
                            )
        print(res.t)
    }
    else if(x$method =="MK"){
        cat(paste("Mann-Kendall Test"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("two-sided homogeinity test"),fill=TRUE)
        cat(paste("H0: S = 0 (no trend)"), fill=TRUE)
        cat(paste("HA: S != 0 (monotonic trend)"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
               cat(paste("Statistics for total series"), fill=TRUE )
        res.t <- data.frame(S = x$Stot,
                            varS = x$Varianz,
                            Z = round(x$Z,1),
                            tau = round(x$tautot,3),
                            pvalue = format.pval(x$pvalue)
                            )
        print(res.t)

    }
#    else if(x$method =="Partial"){
#        cat(paste("Bi-variate partial Mann-Kendall Test"), fill=TRUE)
#        cat(paste(" "), fill=TRUE)
#        cat(paste("two-sided homogeinity test"),fill=TRUE)
#        cat(paste("H0: no trend"), fill=TRUE)
#        cat(paste("HA: monotonic trend"), fill=TRUE)
#        cat(paste(" "), fill=TRUE)
#        cat(paste("Call: "),fill=TRUE)
#        print(x$call)
#        cat(paste(" "), fill=TRUE)
#        cat(paste("Statistics for Variable"), fill=TRUE )
#        res.t <- data.frame(S = x$Stot,
#                            Z = round(x$Z,1),
 #                           tau = round(x$tautot,3),
#                            pvalue = format.pval(x$pvalue)
#                            )
#        print(res.t)
#        cat(paste(" "), fill=TRUE)
#        cat(paste("Correlation matrix: "),fill=TRUE)
#        a <- x$call[[2]]
#        dimnames(x$Correl)[[1]] <- c(a[[2]], a[[3]])
#        dimnames(x$Correl)[[2]] <- c(a[[2]], a[[3]])
#        print(x$Correl)
#    }
    else if(x$method=="SSLP"){
        cat(paste(" "), fill=TRUE)
        cat(paste("Sen's slope and intercept"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("slope: ", round(x$b.sen, 4)), fill=TRUE)
        cat(paste(x$conf.int, "percent confidence intervall for slope"),fill=TRUE)
        cat(paste(round(x$b.sen.up, 4), round(x$b.sen.lo, 4)), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("intercept:", round(x$intercept, 4)), fill=TRUE)
        cat(paste("nr. of observations:", x$nobs), fill=TRUE)
    }
    else if(x$method=="SeaSLP"){
        cat(paste(" "), fill=TRUE)
        cat(paste("Seasonal Sen's slope and intercept"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("slope: ", round(x$b.sen, 4)), fill=TRUE)
        cat(paste("intercept:", round(x$intercept, 4)), fill=TRUE)
        cat(paste("nr. of observations:", x$nobs), fill=TRUE)
    }
}

