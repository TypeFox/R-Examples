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

#' Design: Gaussian responses (known heteroscedasticity)
#'
#' Method for generating a sampling design for data generation
#' following a random effects meta regression model with
#' known heteroscedasticity.
#'
#' Generates a sampling design for the heterogeneity 'h' and
#' a heteroscedasticity 'd1', ..., 'dk'.
#'
#' Points in the design are selected via a maxi-min hypercube sampling
#' using the 'lhs' package in a predefined parameter cube.
#'
#' @param n resolution of the heterogeneity and heteroscedasticity
#' parameters, i.e. the number of of different (heterogeneity,
#' heteroscedasticity) pairs in the design.
#' @param h_bounds bounds of the heterogeneity.
#' @param d_bounds bounds of the heteroscedasticity.
#' @param x design matrix.
#' @return
#' Function returns a data frame.  Each line of this data frame
#' can be an input to the function 'rY' which is used to sample
#' data from such a design.
#' @examples
#' dY <- designY(n=15L, h_bounds=c(0,1), d_bounds=c(0.01,2),
#' x=cbind(1,1:7))
#'
#' if(!all(dim(dY) == c(15,dim(cbind(1,1:7))[1]+1))) {
#'   stop("Wrong dimension")
#' }
#' @export
designY <- function(  n # resolution
                    , h_bounds # bounds of the heterogeneity
                    , d_bounds # bounds of the heteroscedasticity
                    , x # design matrix
                    ) {
    checkArg(n, "integer", len=1, lower=1, na.ok=FALSE)
    checkArg(h_bounds, "numeric", len=2, lower=0, na.ok=FALSE)
    checkArg(d_bounds, "numeric", len=2, lower=0, na.ok=FALSE)
    ps <- makeParamSet(  makeNumericParam("h" , lower=h_bounds[1],
                                          upper=h_bounds[2])
                       , makeNumericVectorParam(  "d"
                                                , len=dim(x)[1]
                                                , lower=d_bounds[1]
                                                , upper=d_bounds[2]))
    return(generateDesign(n=n, ps, fun=maximinLHS))
}

#' Design: Gaussian responses (unknown heteroscedasticity)
#'
#' Method for generating a sampling design for data generation
#' following a random effects meta regression model with
#' unknown heteroscedasticity.
#'
#' Generates a sampling design for the heterogeneity 'h',
#' heteroscedasticity 'd1', ..., 'dk', and study sizes 's1', ..., 'sk'.
#' This design can be used for testing methods that adjust for
#' uncertainty in the heteroscedasticity estimates by additionally
#' considering the size of the respected studies.
#'
#' Points in the design are selected via a maxi-min hypercube sampling
#' using the 'lhs' package in a predefined parameter cube.
#'
#' @param n resolution of the heterogeneity and heteroscedasticity
#' parameters, i.e., the number of of different (heterogeneity,
#' heteroscedasticity, sizes) tuple in the design.
#' @param h_bounds bounds of the heterogeneity.
#' @param d_bounds bounds of the heteroscedasticity.
#' @param s_bounds bounds of the study sizes.
#' @param x design matrix.
#' @return
#' Function returns a data frame.  Each line of this data frame
#' can be an input to the function 'rD' which is used to sample
#' data from such a design.
#' @examples
#' dD <- designD(n=15L, h_bounds=c(0,1), d_bounds=c(0.01,2),
#'   s_bounds=c(200L,2000L), x=cbind(1,1:7))
#'
#' if(!all(dim(dD) == c(15,2*dim(cbind(1,1:7))[1]+1))) {
#'   stop("Wrong dimension")
#' }
#' @export
designD <- function(  n # resolution
                    , h_bounds # bounds of the heterogeneity
                    , d_bounds # bounds of the heteroscedasticity
                    , s_bounds # bounds of the study sizes
                    , x # design matrix
                    ) {
    checkArg(n, "integer", len=1, lower=1, na.ok=FALSE)
    checkArg(h_bounds, "numeric", len=2, lower=0, na.ok=FALSE)
    checkArg(d_bounds, "numeric", len=2, lower=0, na.ok=FALSE)
    checkArg(s_bounds, "integer", len=2, lower=1, na.ok=FALSE)
    ps <- makeParamSet(  makeNumericParam("h" , lower=h_bounds[1],
                                          upper=h_bounds[2])
                       , makeNumericVectorParam(  "d"
                                                , len=dim(x)[1]
                                                , lower=d_bounds[1]
                                                , upper=d_bounds[2]))
    psSizes  <- makeParamSet(makeIntegerVectorParam(  "s"
                                                , len=dim(x)[1]
                                                , lower=s_bounds[1]
                                                , upper=s_bounds[2]))
    return(cbind(  generateDesign(n=n, ps, fun=maximinLHS)
                 , generateDesign(n=n, psSizes, fun=maximinLHS)))
}

#' Design: Binomial responses
#'
#' Method for generating a sampling design for data generation
#' following a binomial-Gaussian model.
#'
#' Generates a sampling design for the heterogeneity 'h', balancing
#' factors 'a1', ..., 'ak' of group assignments, and study sizes 's1',
#' ..., 'sk'.  This design can be used for testing methods for inference
#' for the random effects meta regression model since the logarithm of
#' relative risks of each study is approximately Gaussian distributed.
#' One may use methods that adjust for uncertainty in the
#' heteroscedasticity estimates by additionally considering the size of
#' the respected studies.
#'
#' Points in the design are selected via a maxi-min hypercube sampling
#' using the 'lhs' package in a predefined parameter cube.
#'
#' @param n resolution of the heterogeneity.  n is the number of
#' of different heterogeneity parameters in the design.
#' @param h_bounds bounds of the heterogeneity.
#' @param a_bounds bounds of the balancing factor of group assignments.
#' @param s_bounds bounds of the study sizes.
#' @param r fixed risk in the control.
#' @param x design matrix.
#' @return
#' Function returns a data frame.  Each line of this data frame
#' can be an input to the function 'rB' which is used to sample
#' data from such a design.
#' @examples
#' dB <- designB(n=15L, h_bounds=c(0,1), a_bounds=c(-.3,3),
#'   s_bounds=c(200L,2000L), r=0.03, x=cbind(1,1:5))
#'
#' if(!all(dim(dB) == c(15,2*dim(cbind(1,1:5))[1]+2))) {
#'   stop("Wrong dimension")
#' }
#' @export
designB <- function(  n
                    , h_bounds
                    , a_bounds
                    , s_bounds
                    , r
                    , x
                    ) {
    checkArg(n, "integer", len=1, lower=1, na.ok=FALSE)
    checkArg(h_bounds, "numeric", len=2, lower=0, na.ok=FALSE)
    checkArg(a_bounds, "numeric", len=2, na.ok=FALSE)
    checkArg(s_bounds, "integer", len=2, lower=1, na.ok=FALSE)
    checkArg(r, "numeric", len=1, lower=.0001, upper=.9999)
    ps <- makeParamSet(  makeNumericParam("h"
                                          , lower=h_bounds[1]
                                          , upper=h_bounds[2])
                       , makeNumericVectorParam(  "a"
                                                , len=dim(x)[1]
                                                , lower=a_bounds[1]
                                                , upper=a_bounds[2]))
    psSizes  <- makeParamSet(makeIntegerVectorParam(  "s"
                                                , len=dim(x)[1]
                                                , lower=s_bounds[1]
                                                , upper=s_bounds[2]))
    return(cbind(  generateDesign(n=n, ps, fun=maximinLHS)
                 , generateDesign(n=n, psSizes, fun=maximinLHS)
                 , data.frame(r=r)))
}

#' Data generation: Sampling data of clinical trials
#'
#' A random draw of a hierarchical binomial Gaussian model.
#'
#' It is always assumed that at least one response in a study has
#' happend, i.e., a response of 0 in a treatment or control group is
#' rounded up to 1.  Note that this may lead to an overestimation of
#' small risks.  If possible, make sure your sample sizes are large
#' enough to compensate for this effect.
#'
#' You may work around this by increasing study sizes.
#' @param h heterogeneity.
#' @param s study sizes.
#' @param a balance of group assignments.
#' @param r fixed risk in the control group.
#' @param x design matrix.
#' @param b regression coefficients.
#' @return A list containing the risk and a data frame with the studies.
#' @examples
#' h_test <- .03
#' s_test <- rep(2000, 13)
#' a_test <- rep(.3, 13)
#' x_test <- cbind(1,1:13)
#' b_test <- c(0.02, 0.03)
#' dat <- rBinomGauss(h=h_test, s=s_test, a=a_test, r=0.03 , x=x_test,
#' b=b_test)$study
#'
#' if(!all(dim(dat) == c(dim(x_test)[1], 4))) stop("Wrong dimension")
#' @export
rBinomGauss <- function(  h # heterogeneity
                        , s # study sizes
                        , a # balance of group assignments
                        , r # fixed risk in the control group
                        , x # design matrix
                        , b # regression coefficients
                        ) {
    k  <- length(s)
    y  <- rnorm(k, mean=as.vector(x%*%b), sd=sqrt(h))
    p  <- max(min(r * exp(-y), 0.999), 0.001)
    sizeP <- round(s*(a+1)/2)
    sizeR <- s - sizeP
    n1 <- rbinom(k, sizeP, p)
    m1 <- rbinom(k, sizeR, r)
    dat <- cbind(n1=n1, n2=sizeP-n1, m1=m1, m2=sizeR-m1)
    return(list(  risk=p
                , study=apply(dat, c(1,2), function(tmp) {max(1,tmp)})
                ))
}


#' Data generation: Sampling data of clinical trials
#'
#' Calculates log risk ratios from a study in the right format.
#'
#' @param study Study data of a clinical trial with
#' binomial outcomes.
#' @export
#' @examples
#' h_test <- .03
#' x_test <- cbind(1,1:13)
#' b_test <- c(0.02, 0.03)
#' s_test <- rep(2000, 13)
#' a_test <- rep(.3, 13)
#' rBinomGauss(  h=h_test, s=s_test, a=a_test, r=0.03
#'             , x=x_test, b=b_test)$study -> test
#' yvec(test)
#' dvec(test)
yvec <- function(study) {
    return(log(study[,1]/(study[,1]+study[,2]))
           - log(study[,3]/(study[,3]+study[,4])))
}

#' Data generation: Sampling data of clinical trials
#'
#' Calculates the variance estimate of log risk ratios from a study in
#' the right format.  See the example below for details.
#'
#' @param study Study data of a clinical trial with
#' binomial outcomes.
#' @export
#' @examples
#' h_test <- .03
#' x_test <- cbind(1,1:13)
#' b_test <- c(0.02, 0.03)
#' s_test <- rep(2000, 13)
#' a_test <- rep(.3, 13)
#' rBinomGauss(  h=h_test, s=s_test, a=a_test, r=0.03
#'             , x=x_test, b=b_test)$study -> test
#' yvec(test)
#' dvec(test)
dvec <- function(study) {
    return((1/study[,1]) - (1/(study[,1] + study[,2])) + (1/study[,3]) - (1/(study[,3]+ study[,4])))
}

#' Data generation: Gaussian-Gaussian model
#'
#' Random draws of response vectors y following the distribution
#' of a random effects meta regression model.  Each column is
#' an independent draw.
#'
#' @param n number of draws.
#' @param h heterogeneity.
#' @param d heteroscedasticity.
#' @param x design matrix.
#' @param b regression coefficients.
#' @return A (k,n)-matrix.  Each column is an independent draw.
#' @examples
#' x_test = cbind(1,1:13)
#' h_test = .03
#' d_test = rchisq(13, df=0.02)
#' b_test = c(0.02, 0.03)
#' rY(n=10, h=h_test, d=d_test, x=x_test, b=b_test)
#' @export
rY <- function(  n # number of draws
               , h # heterogeneity
               , d # heteroscedasticity
               , x # design matrix
               , b # regression coefficients
               ) {
    return(matrix(rnorm(  n*length(d)
                        , mean=as.vector(x%*%b)
                        , sd=sqrt(h + d)), ncol=n))
}

#' Data generation: Gaussian-Gaussian model
#'
#' Random draws of heteroscedasticity responses of studies, where each
#' study in a random effects meta regression model follows a Gaussian
#' response. Thus D = (d * X) / (s-1) where X is chi-squared
#' distributed.
#'
#' @param n number of draws.
#' @param d heteroscedasticity.
#' @param s study sizes.
#' @return A (k,n)-matrix.  Each column is an independent draw.
#' @examples
#' d_test = rchisq(13, df=0.02)
#' s_test = rep(100, 13)
#' rD(n=10, d=d_test, s=s_test)
#' @export
rD <- function(  n # number of draws
               , d # heteroscedasticity
               , s # study sizes
               ) {
    return(matrix(d * rchisq(n*length(s), df=(s-1)) / (s-1), ncol=n))
}

#' Data generation: Log-risk-ration of a binomial-Gaussian model
#'
#' Random draws of log risk ratios from a hierarchical binomial
#' Gaussian model.
#'
#' It is always assumed that at least one response in a study has
#' happend, i.e., a response of 0 in a treatment or control group is
#' rounded up to 1.  Note that this may lead to an overestimation of
#' small risks.  If possible, make sure your sample sizes are large
#' enough to compensate for this effect.
#'
#' @param n number of draws.
#' @param h heterogeneity.
#' @param s study sizes.
#' @param a balance of group assignments.
#' @param r fixed risk in the treatment group.
#' @param x design matrix.
#' @param b regression coefficients.
#' @return A (2k,n) matrix.  Each column is an independent draw.
#' @examples
#' h_test <- .03
#' x_test <- cbind(1,1:13)
#' b_test <- c(0.02, 0.03)
#' s_test <- rep(2000, 13)
#' a_test <- rep(.3, 13)
#' rB(n=10, h=h_test, s=s_test, a=a_test, r=.3, x=x_test, b=b_test)
#' @export
rB <- function(  n # number of draws
               , h # heterogeneity
               , s # study sizes
               , a # balance of group assignments
               , r # fixed risk in the treatment group
               , x # design matrix
               , b # regression coefficients
               ) {
    func <- function() {
        study <- rBinomGauss(h=h,s=s,a=a,r=r,x=x,b=b)$study
        return(c(yvec(study=study), dvec(study=study)))
    }

    return(replicate(n, func()))
}
