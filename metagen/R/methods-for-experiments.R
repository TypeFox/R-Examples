# Copyright (C) 2011-2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
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

#' Running a computer experiment
#'
#' Runs a computer experiment that evaluates the performance
#' of different inference methods for the random
#' effects meta regression model with respect to heterogeneity
#' and regression coefficients.
#'
#' @param n number of draws.
#' @param h heterogeneity.
#' @param d heteroscedasticity.
#' @param x design matrix.
#' @param b regression coefficients.
#' @param sgnf significance levels.
#' @param piv_draws privotal draws.
#' @return Data frame of accumulated performance results.
#' @examples
#' h_test <- 0.03
#' x_test <- cbind(1,1:7)
#' b_test <- c(.5, .25)
#' sgnf_test <- c(0.025, 0.01)
#'
#' set.seed(5133568) # for reproducibility
#' d_test <- rchisq(7, df=0.02)
#'
#' # In an actual computer experiment, use 'piv_draws=1000' instead!!
#' experimentY(n=5, h=h_test, d=d_test, x=x_test, b=b_test,
#'   sgnf=sgnf_test, piv_draws=50)
#' @export
experimentY <- function(  n # number of draws
                        , h # heterogeneity
                        , d # heteroscedasticity
                        , x # design matrix
                        , b # regression coefficients
                        , sgnf # significance levels
                        , piv_draws # privotal draws
                        ) {
    samples <- as.data.frame(rY(n=n, h=h, d=d, x=x, b=b))
    accum_results <- function(accum, y) {
        analysis <- metagen(y=y, d=d, x=x, sgnf=sgnf, n=piv_draws,
                            adjusted=F)
        return(Map(rbind, accum, analysis))
    }
    return(Reduce(accum_results, samples, init=metagenEmpty()))
}

#' Running a computer experiment
#'
#' Runs a computer experiment that evaluates the performance
#' of different inference methods for the random
#' effects meta regression model with respect to heterogeneity
#' and regression coefficients.
#'
#' This also includes methods adjusting for uncertainty in
#' the heteroscedasticity vector.  In particular, the study sizes
#' need to be known, here.
#'
#' @param n number of draws.
#' @param h heterogeneity.
#' @param d heteroscedasticity.
#' @param s vector study sizes.
#' @param x design matrix.
#' @param b regression coefficients.
#' @param sgnf significance levels.
#' @param piv_draws privotal draws.
#' @return Data frame of accumulated performance measures.
#' @examples
#' h_test <- 0.03
#' x_test <- cbind(1,1:7)
#' b_test <- c(.5, .25)
#' sgnf_test <- c(0.025, 0.01)
#'
#' set.seed(5133568) # for reproducibility
#' d_test <- rchisq(7, df=0.02)
#' s_test <- runif(7, min=200, max=2000)
#'
#' # In an actual computer experiment, use 'piv_draws=1000' instead!!
#' experimentD(n=5, h=h_test, d=d_test, s=s_test, x=x_test, b=b_test,
#'   sgnf=sgnf_test, piv_draws=50)
#' @export
experimentD <- function(  n # number of draws
                        , h # heterogeneity
                        , d # heteroscedasticity
                        , s # vector study sizes
                        , x # design matrix
                        , b # regression coefficients
                        , sgnf # significance levels
                        , piv_draws # privotal draws
                        ) {
    k <- dim(x)[1]
    samples_y <- rY(n=n, h=h, d=d, x=x, b=b)
    samples_d <- rD(n=n, d=d, s=s)
    samples <- as.data.frame(rbind(samples_y, samples_d))

    accum_results <- function(accum, smpl) {
        smpl_y <- smpl[ (1:k)]
        smpl_d <- smpl[-(1:k)]
        analysis <- metagen(y=smpl_y, d=smpl_d, x=x, sgnf=sgnf,
                            n=piv_draws, s=s, adjusted=T)
        return(Map(rbind, accum, analysis))
    }

    return(Reduce(accum_results, samples, init=metagenEmpty()))
}

### Adding performance measurements to accumulated
### results of a computer experiment running multiple
### analysis of different simulated data following
### a random effects meta regression model.

#' Running a computer experiment: Adding performance measures
#'
#' Adding performance measurements to accumulated
#' results of a computer experiment running multiple
#' analysis of different simulated data following
#' a random effects meta regression model.
#'
#' Adds performance measurements to interval estimates of the
#' regression coefficients.
#' @param accum_int accumulated interval estimates.  At least the
#' following columns need to be present: lower and upper and parameter.
#' @param true true parameter.
#' @examples
#' # For an example, see the 'performance' function.
#' @export
performanceConfR <- function (  accum_int
                              , true
                              ) {
    property <- function (cnf) {
        width <- mean(cnf$upper - cnf$lower, na.rm=T)
        cover <- mean((cnf$lower <= true[cnf$parameter]) &
                        (cnf$upper >= true[cnf$parameter]), na.rm=T)
        return(data.frame(width=width, coverage=cover))
    }
    return(ddply(  .data=accum_int
                 , .variables=grep(  "[^(^lower$)(^upper$)]"
                                   , names(accum_int)
                                   , value=TRUE)
                 , .fun=property))
}

#' Running a computer experiment: Adding performance measures
#'
#' Adding performance measurements to accumulated
#' results of a computer experiment running multiple
#' analysis of different simulated data following
#' a random effects meta regression model.
#'
#' Adds performance measurements to interval estimates of the
#' heterogeneity.
#' @param accum_int accumulated interval estimates.  At least the
#' following columns need to be present: lower and upper.
#' @param true true parameter.
#' @examples
#' # For an example, see the 'performance' function.
#' @export
performanceConfH <- function (  accum_int
                              , true
                              ) {
    property <- function (cnf) {
        width <- mean(cnf$upper - cnf$lower, na.rm=T)
        cover <- mean((cnf$lower <= true) & (cnf$upper >= true),
                      na.rm=T)
        return(data.frame(width=width, coverage=cover))
    }
    return(ddply(  .data=accum_int
                 , .variables=grep(  "[^(^lower$)(^upper$)]"
                                   , names(accum_int)
                                   , value=TRUE)
                 , .fun=property))
}

#' Running a computer experiment: Adding performance measures
#'
#' Adding performance measurements to accumulated
#' results of a computer experiment running multiple
#' analysis of different simulated data following
#' a random effects meta regression model.
#'
#' Adds performance measurements to point estimates of the
#' regression coefficients.
#' @param point accumulated point estimates.
#' @param b true parameter.
#' @examples
#' # For an example, see the 'performance' function.
#' @export
performancePointR <- function (  point # accumulated point estimates.
                               , b     # true coefficients.
                               ) {
    ## Adds performance measurements to point estimates of the
    ## regression coefficients.
    ## Unit test: mse = bias^2 + variance
    property <- function (estimate) {
        difference <- estimate$value - b[estimate$parameter]
        mse      <- mean(difference^2, na.rm=T)
        variance <- var(estimate$value, na.rm=T)
        bias     <- mean(difference, na.rm=T)
        return(data.frame(mse=mse,variance=variance,bias=bias))
    }

    return(ddply(  point
                 , grep("[^(^value$)]", names(point), value=TRUE)
                 , property))
}

#' Running a computer experiment: Adding performance measures
#'
#' Adding performance measurements to accumulated
#' results of a computer experiment running multiple
#' analysis of different simulated data following
#' a random effects meta regression model.
#'
#' Adds performance measurements to point estimates of the
#' heterogeneity.
#' @param point accumulated point estimates.
#' @param h true parameter.
#' @examples
#' # For an example, see the 'performance' function.
#' @export
performancePointH <- function (  point # accumulated point estimates.
                               , h     # true heterogeneity.
                               ) {
    ## Adds performance measurements to point estimates of the
    ## heterogeneity.
    ## Unit test: mse = bias^2 + variance
    property <- function (estimate) {
        difference <- estimate$h - h
        mse      <- mean(difference^2, na.rm=T)
        variance <- var(estimate$h, na.rm=T)
        bias     <- mean(difference, na.rm=T)
        return(data.frame(mse=mse,variance=variance,bias=bias))
    }

    return(ddply(  point
                 , grep("[^(^h$)]", names(point), value=TRUE)
                 , property))
}

#' Running a computer experiment
#'
#' Adding performance measures to the results
#'
#' Calculating performance measurements from a computer experiment.

#' @param results Needs to be of the same type as, for example, the
#' return value of the computer experiments 'experimentY',
#' 'experimentD'.
#' @param b true regression coefficients.
#' @param h true heterogeneity.
#' @return Data frame containing performance measurements of inference
#' methods based on the results of the computer experiment given by
#' 'results'.
#' @examples
#' h_test <- 0.03
#' x_test <- cbind(1,1:7)
#' b_test <- c(.5, .25)
#' sgnf_test <- c(0.025, 0.01)
#'
#' set.seed(5133568) # for reproducibility
#' d_test <- rchisq(7, df=0.02)
#'
#' # In an actual computer experiment, use 'piv_draws=1000' instead!!
#' eY <- experimentY(n=5, h=h_test, d=d_test, x=x_test, b=b_test,
#'   sgnf=sgnf_test, piv_draws=50)
#'
#' performance(results=eY, b=b_test, h=h_test)
#' @export
performance <- function (  results # results of a computer experiment
                         , b       # true regression coefficients
                         , h       # true heterogeneity
                         ) {
    confr <- performanceConfR(results$confr, true=b)
    confh <- performanceConfH(results$confh, true=h)
    pointr <- performancePointR(results$pointr, b=b)
    pointh <- performancePointH(results$pointh, h=h)
    return(list(pointh=pointh, confh=confh, pointr=pointr, confr=confr))
}
