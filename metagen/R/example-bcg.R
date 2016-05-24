# Copyright (C) 2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
#
#     This program is free software: you can redistribute it and/or
#     modify it under the terms of the GNU General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see
#     <http://www.gnu.org/licenses/>.

#' Example: Setting up the BCG-data set
#'
#' Exemplary data set of 14 clinical trials evaluating BCG vaccine
#' efficacy.
#'
#' Reads in the BCG vaccine efficacy data from the metafor package and
#' adds some statistics to the data such as the log-relative risk, study
#' size, measurements of balance, confidence intervals of the responses,
#' and the like.
#'
#' @param sgnf significance level of the confidence intervals for the
#' relative risks.
#' @return Returns a data set of 13 clinical trials which evaluated the
#' efficacy of the BCG vaccine.  The data set is an exact copy of the
#' data set found in the dat.bcg data frame provided by the metafor
#' package.
#' @examples
#' bcgVaccineData()
#' @export
bcgVaccineData <- function (sgnf=0.05) {
    data("dat.bcg", package = "metafor", envir = environment())
    # RR=log relative risk
    dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
                  data=dat.bcg, append=TRUE)

    tmp_size_t  = dat$tpos + dat$tneg
    tmp_size_c  = dat$cpos + dat$cneg
    tmp_size    = tmp_size_t + tmp_size_c
    tmp_balance = (tmp_size_t - tmp_size_c) / tmp_size

    bcg_dat <- data.frame(
      reference= paste(gsub("et al", "et al.", dat$author), dat$year, sep=", ")
    , x=         dat$ablat
    , size=      tmp_size
    , balance=   tmp_balance
    , logrisk=   unname(dat$yi) #
    , sdiv=      dat$vi
    , rlower=    dat$yi - qnorm(1-(sgnf/2)) * sqrt(dat$vi)
    , rupper=    dat$yi + qnorm(1-(sgnf/2)) * sqrt(dat$vi)
    , allocation= factor(dat$alloc)
    )

    ## Order the date frame by absolute latitude
    bcg_dat <- with(bcg_dat, bcg_dat[order(x),])
    bcg_dat$reference <- with(bcg_dat,
        factor(order(x), labels=paste(reference[order(x)]), ordered=T))
    rownames(bcg_dat) <- NULL
    return(bcg_dat)
}

#' Example: Plotting study sizes
#'
#' @param dat data frame of study responses of binomial type.
#' @return
#' An object created by ggplot2.
#' @examples
#' bcg <- bcgVaccineData()
#' plotStudySizes(bcg)
#' @export
plotStudySizes <- function(dat) {
    return(ggplot(dat, aes(reference, weight=size))
           + geom_bar()
           + scale_x_discrete("")
           + scale_y_continuous(expression("Total number of subjects"))
           + coord_flip()
           + theme(axis.text = element_text(colour = "black"))
           )
}

#' Example: Plotting study unbalances in group assignments
#'
#' @param dat data frame of study responses of binomial type.
#' @return
#' An object created by ggplot2.
#' @export
plotStudyUnbalance <- function(dat) {
    return( ggplot(dat, aes(reference, weight=balance))
           + geom_bar()
           + scale_x_discrete("")
           + scale_y_continuous(expression(frac(v-neg(v), v+neg(v))))
           + coord_flip()
           + theme(axis.text = element_text(colour = "black"))
           )
}

#' Example: Plotting a forest plot of a data frame
#'
#' @param dat data frame of study responses of binomial type.
#' @return
#' An object created by ggplot2.
#' @examples
#' bcg <- bcgVaccineData()
#' plotStudyForest(bcg)
#' @export
plotStudyForest <- function(dat) {
    return( ggplot(dat, aes(reference, y=logrisk, ymin=rlower,
                            ymax=rupper))
           + geom_pointrange(shape=15, size=1)
           + geom_hline(yintercept=0, size=.3)
           + scale_x_discrete("")
           + scale_y_continuous("Logarithm of relative risk")
           + coord_flip()
           + theme(axis.text = element_text(colour = "black"))
           )
}

#' Example: Plotting the q- and p-function from the dissertation
#'
#' @param y a vector of responses.
#' @param d a vector of heteroscedasticity.
#' @param x a design matrix.
#' @param n number of points to interpolate along.
#' @return
#' A list of objects created by ggplot2.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#' p <- plotStudyQfuncPfunc(y=bcg_y, d=bcg_d, x=bcg_x, n=500)
#' p[1] # plot of the q-function
#' p[2] # plot of the p-funciton
#' @export
plotStudyQfuncPfunc <- function (  y # a vector of responses.
                                 , d # a vector of heteroscedasticity.
                                 , x # a design matrix.
                                 , n # number of points to interpolate
                                     # along.
                                 ) {
    qfu <- qfunc(y,d,x)
    pfu <- pfunc(y,d,x)

    findk <- 0
    while (qfu(exp(findk)) > 0.1*qfu(0)) {findk=findk + 1}

    plotQfunc <- (  qplot(c(0,exp(findk)), c(0,0), geom="blank")
                  + stat_function(fun = qfu, colour="brown", n=n)
                  + geom_hline(yintercept=0)
                  + geom_vline(xintercept=0)
                  + scale_x_continuous(expression(tau))
                  + scale_y_continuous(expression(q[delta](tau)))
                  + theme(axis.text = element_text(colour = "black"))
                  )

    plotPfunc1 <- ( qplot(c(0,exp(findk)), c(0,0), geom="blank")
                   + stat_function(fun = qfu, colour="brown", n=n)
                   + geom_hline(yintercept=0)
                   + geom_vline(xintercept=0)
                   + scale_x_continuous(expression(p[delta](eta)))
                   + scale_y_continuous(expression(eta))
                   + coord_flip()
                   + theme(axis.text = element_text(colour = "black"))
                   )

    plotPfunc2 <- ( qplot(c(0,qfu(0)), c(0,0), geom="blank")
                   + stat_function(fun = pfu, colour="brown", n=n)
                   + geom_hline(yintercept=0)
                   + geom_vline(xintercept=0)
                   + scale_x_continuous(expression(eta))
                   + scale_y_continuous(expression(p[delta](eta)))
                   + theme(axis.text = element_text(colour = "black"))
                   )

    return(list(plotQ=plotQfunc, plotP=plotPfunc1,
                plotP_debug=plotPfunc2))
}

#' Example: Plotting interval estimates
#'
#' Plots a graphical representation of interval estimates
#' in the data frame 'cnf' by type of method used for
#' the estimation.
#' @param cnf data frame of interval estimates
#' @return
#' An object created by ggplot2.
#' @export
plotIntervalEstimates <- function (cnf) {
    return(ggplot(  cnf
                  , aes(type, 0, ymin = lower, ymax=upper))
    + geom_linerange(size=5)
    + scale_x_discrete(name="")
    + geom_hline(yintercept=0, size=.3)
    + coord_flip()
    + facet_grid(confidence ~ .)
    + ylab("Density")
    + theme(axis.text = element_text(colour = "black"))
    )
}
