summary.trend.test <-
function (object, ...){
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
##    This is a summary function.
##
    if (object$method =="SMK"){
        cat(paste("Seasonal Mann-Kendall Test without correlation"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("two-sided homogeinity test"),fill=TRUE)
        cat(paste("H0: S = 0 (no trend)"), fill=TRUE)
        cat(paste("HA: S != 0 (monotonic trend)"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("Statistics for individual seasons"), fill=TRUE)
        res.S <- data.frame(S = object$Sg,
                            varS = round(object$varSg,1),
                            Z = round(object$Zg,1),
                            tau = round(object$taug,3),
                            pvalue = format.pval(object$pvalg)
                            )

        print(res.S)
        cat(paste(" "), fill=TRUE)
        cat(paste("Statistics for total series"), fill=TRUE )
        res.t <- data.frame(S = object$Stot,
                            varS = object$Varianz,
                            Z = round(object$Z,1),
                            tau = round(object$tautot,3),
                            pvalue = format.pval(object$pvalue)
                            )
        print(res.t)
    }
    else if(object$method =="CSMK"){
        cat(paste("Seasonal Mann-Kendall Test with correlations"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("two-sided homogeinity test"),fill=TRUE)
        cat(paste("H0: S = 0 (no trend, seasons are correlated)"), fill=TRUE)
        cat(paste("HA: S != 0 (monotonic trend)"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
               cat(paste("Statistics for total series"), fill=TRUE )
        res.t <- data.frame(S = object$Stot,
                            varS = object$Varianz,
                            Z = round(object$Z,1),
 #                           tau = round(object$tautot,3),
                            pvalue = format.pval(object$pvalue)
                            )
        print(res.t)
    }
    else if(object$method =="MK"){
        cat(paste("Mann-Kendall Test"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("two-sided homogeinity test"),fill=TRUE)
        cat(paste("H0: S = 0 (no trend)"), fill=TRUE)
        cat(paste("HA: S != 0 (monotonic trend)"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
               cat(paste("Statistics for total series"), fill=TRUE )
        res.t <- data.frame(S = object$Stot,
                            varS = object$Varianz,
                            Z = round(object$Z,1),
                            tau = round(object$tautot,3),
                            pvalue = format.pval(object$pvalue)
                            )
        print(res.t)

    }
#    else if(object$method =="Partial"){
#        cat(paste("Bi-variate partial Mann-Kendall Test"), fill=TRUE)
#        cat(paste(" "), fill=TRUE)
#        cat(paste("two-sided homogeinity test"),fill=TRUE)
#        cat(paste("H0: no trend"), fill=TRUE)
#        cat(paste("HA: monotonic trend"), fill=TRUE)
#        cat(paste(" "), fill=TRUE)
#        cat(paste("Call: "),fill=TRUE)
#        print(object$call)
#        cat(paste(" "), fill=TRUE)
#        cat(paste("Statistics for Variable"), fill=TRUE )
#        res.t <- data.frame(S = object$Stot,
#                            Z = round(object$Z,1),
 #                           tau = round(object$tautot,3),
#                            pvalue = format.pval(object$pvalue)
#                            )
#        print(res.t)
#        cat(paste(" "), fill=TRUE)
#        cat(paste("Correlation matrix: "),fill=TRUE)
#        a <- object$call[[2]]
#        dimnames(object$Correl)[[1]] <- c(a[[2]], a[[3]])
#        dimnames(object$Correl)[[2]] <- c(a[[2]], a[[3]])
#        print(object$Correl)
#    }
    else if(object$method=="SSLP"){
        cat(paste(" "), fill=TRUE)
        cat(paste("Sen's slope and intercept"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("slope: ", round(object$b.sen, 4)), fill=TRUE)
        cat(paste(object$conf.int, "percent confidence intervall for slope"),fill=TRUE)
        cat(paste(round(object$b.sen.up, 4), round(object$b.sen.lo, 4)), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("intercept:", round(object$intercept, 4)), fill=TRUE)
        cat(paste("nr. of observations:", object$nobs), fill=TRUE)
    }
    else if(object$method=="SeaSLP"){
        cat(paste(" "), fill=TRUE)
        cat(paste("Seasonal Sen's slope and intercept"), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste(" "), fill=TRUE)
        cat(paste("slope: ", round(object$b.sen, 4)), fill=TRUE)
        cat(paste("intercept:", round(object$intercept, 4)), fill=TRUE)
        cat(paste("nr. of observations:", object$nobs), fill=TRUE)
    }
}

