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

## The following two functions should always
## return the same value as 'qfunc'.

qfunc1 <- function (  y   # study responses
                    , d   # heteroscedasticity
                    , x   # design matrix
                    ) {
    ## Returns the q-function. (slow)
    qfuncuno <- function (tau) {
        o <- ginv(diag(tau + d))
        e <- (  diag(1,dim(x)[1]) - x %*%
              ginv(t(x) %*% o %*% x) %*%
              t(x) %*% o)
        return(as.vector( y %*% o %*% e %*% y ))
    }
    return(function (tauvec) {return(sapply(tauvec, qfuncuno))})
}

qfunc2 <- function (  y   # study responses
                    , d   # heteroscedasticity
                    , x   # design matrix
                    ) {
    ## Returns the q-function. (faster)
    e     <- diag(1,dim(x)[1]) - x %*% ginv(t(x) %*% x) %*% t(x)
    eig   <- eigen(e)
    K     <- eig$vectors[,which(eig$values > .5)]
    qfuncuno <- function(tau) {
        oinv <- diag(tau + d)
        tmp  <- t(K) %*% y
        return(as.vector(
                 t(tmp) %*% ginv(t(K) %*% diag(tau+d) %*% K) %*% tmp))
    }
    return(function (tauvec) {return(sapply(tauvec, qfuncuno))})
}

#' The q_delta(tau) function.
#'
#' Returns the q-function.
#'
#' @param y study responses.
#' @param d heteroscedasticity.
#' @param x design matrix.
#' @return A vector valued function.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#' qfunc(y=bcg_y, d=bcg_d, x=bcg_x)
#' @export
qfunc <- function (  y # study responses
                   , d # heteroscedasticity
                   , x # design matrix
                   ) {
    SVD <- svd(t(x),nv=dim(x)[1])
    K   <- (SVD$v)[,-(1:length(SVD$d))]
    tmp <- y%*%K
    qfuncuno <- function (h) {
        return(as.vector(
                 tmp %*% ginv(t(K) %*% diag(h+d) %*% K) %*% t(tmp)))
    }
    return(function (vec) {return(sapply(vec, qfuncuno))})
}

#' The p_delta(eta) function.
#'
#' Returns the p-function.
#'
#' @param y study responses.
#' @param d heteroscedasticity.
#' @param x design matrix.
#' @return A vector valued function.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#' pfunc(y=bcg_y, d=bcg_d, x=bcg_x)
#'
#' # Calculating the Mandel-Paule estimate:
#' pfunc(y=bcg_y, d=bcg_d, x=bcg_x)(dim(bcg_x)[1] - dim(bcg_x)[2])
#' @export
pfunc <- function ( y   # study responses
                  , d   # heteroscedasticity
                  , x   # design matrix
                  ) {
    qfun  <- qfunc(y,d,x)
    qfu0  <- qfun(0)
    tmpf  <- function(tmp, eta) {return(qfun(tmp) - eta)}
    findk <- 0
    pfun  <- function(eta) {
        findk <- 0
        while (tmpf(exp(findk), eta) > 0) {findk = findk+1}
        return(uniroot(  tmpf
                       , c(0,exp(findk))
                       , f.lower=qfu0, eta=eta)$root)
    }
    pfuncuno <- function(eta) {
        tmp=Inf
        if (eta >= qfu0) {tmp=0} else if (eta > 0) {tmp=pfun(eta)}
        return(tmp)
    }
    return(function(etavec) {return(sapply(etavec, pfuncuno))})
}

#' Regression coefficients: formulaL
#'
#' Calculate pivotal quantities for the regression coefficients
#' using the method: formulaL form the dissertation.
#'
#' Algorithm for calculating a single generalised pivotal quantity
#' for the regression coefficients for given generalised pivotal
#' quantities for the heterogeneity using the univariate version
#' of the pivotal formula.
#'
#' @param y k-vector of responses.
#' @param d k-vector of heteroscedasticity.
#' @param h scalar of heterogeneity.
#' @param g p-vector of some p-variate Gaussian draw.
#' @param x design k-p-matrix.
#' @return A p-vector.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#'
#' # When for example using the Mandel-Paule estimate:
#' bcg_h <- pfunc(y=bcg_y, d=bcg_d, x=bcg_x)(dim(bcg_x)[1] -
#'   dim(bcg_x)[2])
#'
#' set.seed(51351) # for reproducibility
#' random_g <- rnorm(dim(bcg_x)[2])
#' formulaL(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=bcg_x)
#'
#' # The function can also be used when planing to perform
#' # a meta regression with no intercept, and only a singel
#' # covariate (i.e. dim(x) = 1).  In this case,
#' # the design matrix can simply be provided by a vector.
#' set.seed(51351) # for reproducibility
#' random_g <- rnorm(1)
#' formulaL(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=bcg$x)
#'
#' # When performing a meta analysis, provide the function
#' # with a vector of 1s.
#' formulaL(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=rep(1,
#'   length(bcg_y)))
#' @export
formulaL <- function ( y # k-vector of responses
                     , d # k-vector of heteroscedasticities
                     , h # scalar of heterogeneity
                     , g # p-vector of some p-variate Gaussian draw
                     , x # design k-p-matrix
                     ) {
    O  <- diag(1/(h+d))
    V  <- t(x) %*% O %*% x
    By <- as.vector(solve(t(x) %*% O %*% x, t(x) %*% O %*% y))
    L  <- By - sqrt(diag(ginv(V))) * g
    return(L)
}

#' Regression coefficients: formulaR
#'
#' Calculate pivotal quantities for the regression coefficients
#' using the method: formulaR form the dissertation.
#'
#' Algorithm for calculating a single generalised pivotal quantity for
#' the regression coefficients for given generalised pivotal quantities
#' for the heterogeneity using the multivariate version of the pivotal
#' formula.
#'
#' @param y k-vector of responses.
#' @param d k-vector of heteroscedasticity.
#' @param h scalar of heterogeneity.
#' @param g p-vector of some p-variate Gaussian draw.
#' @param x design k-p-matrix.
#' @return A p-vector.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#'
#' # When, for example, using the Mandel-Paule estimate:
#' bcg_h <- pfunc(y=bcg_y, d=bcg_d, x=bcg_x)(dim(bcg_x)[1] -
#'   dim(bcg_x)[2])
#'
#' set.seed(51351) # for reproducibility
#' random_g <- rnorm(dim(bcg_x)[2])
#' formulaR(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=bcg_x)
#'
#' # The function can also be used when planing to perform
#' # a meta regression with no intercept, and only a singel
#' # covariate (i.e. dim(x) = 1).  In this case,
#' # the design matrix can simply be provided by a vector.
#' set.seed(51351) # for reproducibility
#' random_g <- rnorm(1)
#' formulaR(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=bcg$x)
#'
#' # When performing a meta analysis, provide the function
#' # with a vector of 1s.
#' formulaR(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=rep(1,
#'   length(bcg_y)))
#' @export
formulaR <- function ( y # k-vector of responses
                     , d # k-vector of heteroscedasticities
                     , h # scalar of heterogeneity
                     , g # p-vector of some p-variate Gaussian draw
                     , x # design k-p-matrix
                     ) {
    ## Algorithm for calculating a single generalised pivotal quantity
    ## for the regression coefficents for given generalised pivotal
    ## quantities for the heterogeneity using the univariate version
    ## of the pivotal formula.
    ##
    ## Return is a p-vector.
    if (is.vector(x)) {x <- as.matrix(x, ncol=1)}
    O <- diag(1/(h+d))
    e <- eigen(t(x) %*% diag(1/(h+d)) %*% x)
    principleRoot <- ( e$vectors %*%
                      diag(sqrt(e$values), ncol=dim(x)[2]) %*%
                      t(e$vectors))
    tranformedY <- as.vector(t(x) %*% diag(1/(h+d)) %*% y)
    tmp <- solve(principleRoot, tranformedY)
    R   <- solve(principleRoot, tmp - g)
    return(R)
}

#' Steams of pivotal quantities of the regression coefficient
#'
#' Algorithm for generating a steam of generalised pivotal quantities
#' for the regression coefficients.  If adjusted=FALSE, then no
#' adjustments are made for the uncertainty in the heteroscedasticity
#' estimates d.  If adjusted=TRUE, then adjustments are performed.  In
#' this case, 's' needs to be provided.
#'
#' @param n      length of stream.
#' @param y      k-vector of responses.
#' @param d      k-vector of heteroscedasticity.
#' @param x      design (k,p)-matrix.
#' @param s      k-vector of study responses.  No need to provide this,
#' when adjusted=FALSE.  Default is NULL.
#' @param method A list.  Used to choose the methods for calculating
#' the pivotal quantities of the regression coefficients.  Default
#' is 'method=list("univariate", "multivariate")'.
#' @param adjusted TRUE or FALSE.  Default is FALSE.
#' @return If method=="univariate" or method=="multivariate", then the
#' return is a (p+1)-n-matrix.  The first row contains pivotal
#' quantities of the heterogeneity, the rest of the rows pivotal
#' quantities of the regression coefficients. Each column is an
#' independent draw.
#'
#' If 'method==list("univariate", "multivariate")', then the return is a
#' (2p+1)-n-matrix.  Of each column, the first element is a pivotal for
#' the heterogeneity, the next 'p' elements is a pivotal vector for the
#' regression coefficients based on "univariate", the last 'p' elements
#' are a pivotal vector for the regression coefficients based on
#' "multivariate"
#' @export
pivotalStream <- function (  n      # length of stream.
                           , y      # k-vector of responses.
                           , d      # k-vector of heteroscedasticities.
                           , x      # design k-p-matrix.
                           , s=NULL # k-vector of study responses
                           , method=list("univariate", "multivariate")
                           , adjusted    # TRUE or FALSE
                           ) {
    if ( ("univariate" %in% method) && !("multivariate" %in% method)) {
        func <- function(  ph # pivotal of the heterogeneity
                         , pd # pivotal of the heteroscedasticity.
                         , rg # a random draw of p-variate Gaussian
                         ) {
            return(formulaL(y, pd, ph, rg, x))
        }
    } else if (!("univariate" %in% method) &&
               ("multivariate" %in% method)) {
        func <- function(  ph # pivotal of the heterogeneity
                         , pd # pivotal of the heteroscedasticity.
                         , rg # a random draw of p-variate Gaussian
                         ) {
            return(formulaR(y, pd, ph, rg, x))
        }
    } else
        func <- function(  ph # pivotal of the heterogeneity
                         , pd # pivotal of the heteroscedasticity.
                         , rg # a random draw of p-variate Gaussian
                         ) {
            return(  c(formulaL(y, pd, ph, rg, x)
                     , formulaR(y, pd, ph, rg, x)))
    }

    if (!adjusted) {
        piv_d <- d
        piv_h <- pfunc(y,d,x)(rchisq(n, dim(x)[1] - dim(x)[2]))
    } else {
        piv_d <- d * (s-1) / matrix(  rchisq(n*dim(x)[1],df=s-1)
                                    , nrow=dim(x)[1])
        piv_h <- apply(piv_d, 2, function (pd) {
                       pfunc(y,pd,x)(rchisq(1,dim(x)[1] - dim(x)[2]))}
        )
    }
    return(rbind(piv_h, mapply(  func
                  , piv_h
                  , as.data.frame(piv_d)
                  , as.data.frame(matrix(  rnorm(n*dim(x)[2])
                                         , nrow=dim(x)[2]))
                  ))
    )
}

#' Inference: Based on generalised inference principles.
#'
#' @param y      k-vector of responses.
#' @param d      k-vector of heteroscedasticities.
#' @param x      design k-p-matrix.
#' @param sgnf   vector of significance levels
#' @param s      k-vector of study responses.  No need to provide this,
#' when 'adjusted==FALSE'.  Default is NULL.
#' @param n      draws from the pivotal distribution.
#' @param method Default is 'list("univariate", "multivariate")'.
#' @param adjusted TRUE or FALSE.  Default is FALSE.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#' sgnf_lev <- c(0.01, 0.025, 0.05, 0.01)
#'
#' set.seed(865287113) # for reproducibility
#'
#' # Runs a standard analysis, use n=1000 in an actual
#' # analysis instead!!
#' g1 <- metagenGeneralised(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=0.025, n=50)
#' g2 <- metagenGeneralised(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=sgnf_lev,
#'   n=50)
#'
#' # Runs the methods based on generalised principles via an
#' # adjustment for the unknown heteroscedasticity.  Use n=1000 in an
#' # actual analysis instead!!
#' bcg_s <- bcg$size
#' g3 <- metagenGeneralised(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=0.025,
#'   s=bcg_s, n=50, adj=TRUE)
#' g4 <- metagenGeneralised(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=sgnf_lev,
#'   s=bcg_s, n=50, adj=TRUE)
#'
#' # The implementation can also handle the case in which
#' # a meta regression is planed with no intercept and only a
#' # single covariate (i.e. dim(x) = 1).  In this case,
#' # the design matrix can simply be provided by a vector.
#' # (This makes no sense in this example and shall only proves
#' # feasibility)
#' g5 <- metagenGeneralised(y=bcg_y, d=bcg_d, x=bcg$x, sgnf=0.025, n=50)
#'
#' # When performing a meta analysis, provide the function
#' # with a vector of 1s.
#' g6 <- metagenGeneralised(y=bcg_y, d=bcg_d, x=rep(1,length(bcg_y)),
#'   sgnf=0.025, n=50)
#'
#' if (!all(names(g1) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(g2) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(g3) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(g4) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(g5) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(g6) == names(metagenEmpty()))) stop("Name clash")
#' @export
metagenGeneralised <- function (  y      # k-vector of responses.
                                , d      # k-vector of heteroscedasticit
                                , x      # design k-p-matrix.
                                , sgnf   # vector of significance levels
                                , s=NULL # k-vector of study responses
                                , n      # draws of pivotal distribution
                                , method=list("univariate", "multivariate")
                                , adjusted=FALSE # TRUE or FALSE
                                ) {
    if (is.vector(x)) {x <- as.matrix(x, ncol=1)}
    pStream <- pivotalStream(n=n, y=y, d=d, s=s, x=x, method=method,
                             adjusted=adjusted)
    quan <- apply(pStream, 1, quantile, probs=c(0.5, sgnf/2, 1-sgnf/2))
    colnames(quan) <- NULL
    rownames(quan) <- NULL
    pointh_values <- quan[1,1]
    pointr_values <- quan[1,-1]
    lower_estimates <- quan[2:(length(sgnf)+1),]
    upper_estimates <- quan[-(1:(length(sgnf)+1)),]
    if (length(sgnf) > 1) {
        bounds_h <- cbind(  as.vector(lower_estimates[,1])
                          , as.vector(upper_estimates[,1]))
        bounds_r <- cbind(  as.vector(lower_estimates[,-1])
                          , as.vector(upper_estimates[,-1]))
    } else {
        bounds_h <- cbind(  as.vector(lower_estimates[1])
                          , as.vector(upper_estimates[1]))
        bounds_r <- cbind(  as.vector(lower_estimates[-1])
                          , as.vector(upper_estimates[-1]))
    }
    colnames(bounds_h) <- c("lower", "upper")
    colnames(bounds_r) <- c("lower", "upper")

    type_str <- if (adjusted) "adjusted" else "unadjusted"
    meth_str <- "generalised"

    confr_names <- data.frame( type=factor(type_str)
        , method=factor(paste(  meth_str
                              , rep(method,
                                    each=(length(sgnf)*(dim(x)[2])))))
        , parameter=factor(rep(1:dim(x)[2], each=length(sgnf)))
        , confidence=1-sgnf)

    confh_names <- data.frame(  type=factor(paste(meth_str, type_str))
                              , confidence=1-sgnf)

    pointr <- data.frame(
        type=factor(paste( meth_str
                          , rep(method, each=(dim(x)[2]))
                          , type_str))
        , parameter=factor(1:dim(x)[2])
        , value=pointr_values)

    pointh <- data.frame(  type=factor(paste(meth_str, type_str))
                         , h=pointh_values)

    return(list(  pointh=pointh
                , confh=cbind(confh_names, bounds_h)
                , pointr=pointr
                , confr=cbind(confr_names, bounds_r)))
}

### Point estimates for the heterogeneity parameter
### and the regression coefficients
###

tryFunc <- function ( func ) {
    ## This is a simple helper function that is used
    ## since some of the iterative algorithms
    ## may not converge in all cases.
    result <- try(func)
    return(if (inherits(result, "try-error")) NA else result)
}

#' Point estimates: For the heterogeneity parameter
#'
#' Returns a list of tau estimates based on different approximative
#' methods.

#' Different point estimates for the heterogeneity parameter are
#' calculated:
#' HD    (Hedges),
#' SL    (DerSimonian-Laird),
#' SJ    (Sidik-Jonkman),
#' MP    (Mandel-Paule),
#' ML    (maximum likelihood),
#' REML  (restricted maximum-likelihood).
#' Since any of these methods may fail to converge,
#' there result may be 'NA' in this case.
#'
#' @param y study responses
#' @param d heteroscedasticity
#' @param x design matrix
#' @return A data frame containing point estimates.  Variables
#' are 'type' and 'h'.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#' hEstimates(y=bcg_y, d=bcg_d, x=bcg_x)
#'
#' # The implementation can also handle the case in which
#' # a meta regression is planed with no intercept and only a
#' # single covariate (i.e. dim(x) = 1).  In this case,
#' # the design matrix can simply be provided by a vector.
#' # (This makes no sense in this example and shall only prove
#' # feasibility)
#' hEstimates(y=bcg_y, d=bcg_d, x=bcg$x)
#'
#' # When performing a meta analysis, provide the function
#' # with a vector of 1s.
#' hEstimates(y=bcg_y, d=bcg_d, x=rep(1, length(bcg_y)))
#' @export
hEstimates <- function (  y    # study responses
                        , d    # heteroscedasticity
                        , x    # design matrix
                        ) {
    if (is.vector(x)) x <- as.matrix(x, ncol=1)
    tauHE   <- data.frame(type="Hedges"
                          , h=tryFunc(rma(y, d, mods=x, method="HE"
                                          , intercept=FALSE)$tau2))
    tauDL   <- data.frame(type="DerSimonian-Laird"
                          , h=tryFunc(rma(y, d, mods=x, method="DL"
                                          , intercept=FALSE)$tau2))
    tauSJ   <- data.frame(type="Sidik-Jonkman"
                          , h=tryFunc(rma(y, d, mods=x, method="SJ"
                                          , intercept=FALSE)$tau2))
    tauML   <- data.frame(type="maximum-likelihood"
                          , h=tryFunc(rma(y, d, mods=x, method="ML"
                                          , intercept=FALSE)$tau2))
    tauREML <- data.frame(type="restricted maximum-likelihood"
                          , h=tryFunc(rma(y, d, mods=x, method="REML"
                                          , intercept=FALSE)$tau2))
    tauMP <- data.frame(type="Mandel-Paule",
                        h=pfunc(y=y,d=d,x=x)(dim(x)[1] - dim(x)[2]))
    return(rbind(tauHE,tauDL,tauSJ,tauMP,tauML,tauREML))
}

coeffEstimates <- function (  y # k-vector of responses
                            , d # k-vector of heteroscedasticities
                            , h # scalar of heterogeneity
                            , x # k-p-matrix (design matrix)
                            ) {
    ## Calculate a point estimate for the regression coefficient
    ## for given point estimates of the variance components 'd' and 'h'.
    O <- diag(1/(h+d))
    V <- t(x) %*% O %*% x
    return(as.vector(solve(t(x) %*% O %*% x, t(x) %*% O %*% y)))
}

#' Point estimates: For the regression coefficients
#'
#' Calculates point estimates for the regression coefficient
#' for given point estimates of the variance components 'd' and a data
#' frame of different estimates of the heterogeneity 'h'.
#'
#' @param y study responses, k-vector of responses.
#' @param d heteroscedasticity, k-vector of heteroscedasticities.
#' @param h_dat Here, 'h_dat' should be a data frame with variables
#' 'type' and 'h'.  Thus, one may use h_dat = hEstimates(y, d, x).
#' @param x design matrix, k-p-matrix.
#' @return A list of estimates for the regression coefficients.
#'
#' Here, 'h_dat' should be a data frame with variables 'type' and 'h',
#' thus, we may use h_dat = hEstimates(y, d, x)
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#' bcg_h <- hEstimates(y=bcg_y, d=bcg_d, x=bcg_x)
#' regressionEstimates(y=bcg_y, d=bcg_d, h_dat=bcg_h, x=bcg_x)
#' @export
regressionEstimates <- function (  y
                                 , d
                                 , h_dat
                                 , x
                                 ) {
    func <- function (r) {
        if (is.na(r$h)) data.frame() else {
            data.frame(  parameter=factor(1:dim(x)[2])
                       , value=coeffEstimates(y=y,d=d,h=r$h,x=x))
        }
    }
    return(ddply(h_dat,"type",func))
}

### Functions for symmetric confidence intervals
### based on standard deviations, point
### estimates and quantile functions
###

#' Interval estimates: Generic function
#'
#' Generic function to produce interval estimates of univariate
#' parameters based on first order limit theory.
#'
#' Function for symmetric confidence intervals based on standard
#' deviations, point estimates, and quantile functions.
#'
#' Can only handle a single significance level!  See
#' 'makeConfInts' for a more flexible solution.
#' @param sgn   one significance level.
#' @param pst   point estimate.
#' @param fct   standard error.
#' @param crt   function for critical value computation.
#' @param name  string: name of the method.
#' @export
makeConfInt <- function(  sgn  # one significance level
                        , pst  # point estimate
                        , fct  # standard error
                        , crt  # function for critical value computation
                        , name # string: name of the method
                        ) {
    ## Returns a confidence interval based on the point
    ## estimate standard error, critical value and the significance
    ## levels given.
    crtval <- crt(sgn/2)
    return(data.frame(  method=name
                      , parameter=factor(1:length(pst))
                      , confidence=1-sgn
                      , lower=pst + fct * crtval
                      , upper=pst - fct * crtval))
}

#' Interval estimates: Generic function
#'
#' Generic function to produce interval estimates of
#' univariate parameters based on first order limit theory.
#'
#' Function for symmetric confidence intervals based on standard
#' deviations, point estimates, and quantile functions.
#'
#' @param sgn   one significance level.
#' @param pst   point estimate.
#' @param fct   standard error.
#' @param crt   function for critical value computation.
#' @param name  string: name of the method.
#' @export
makeConfInts <- function(  sgn  # vector of significance levels
                         , pst  # point estimate
                         , fct  # standard error
                         , crt  # critical value
                         , name # string: name of method
                         ) {
    ## Returns confidence intervals based on the point
    ## estimate standard error, critical value and the significance
    ## levels given.
    return(ldply(sgn, makeConfInt, pst, fct, crt, name))
}

### Calculating stochastic approximations
###

intervalEstimates_OneSgnf <- function (  y    # study responses
                                       , d    # heteroscedasticity
                                       , h    # point estimate of tau
                                       , x    # design matrix
                                       , sgnf # significance levels
                                       ) {
    ## Confidence intervals based on first order limit theory.
    ## TODO: I could also put the point estimates in the return.
    By    <- coeffEstimates(y=y,d=d,h=h,x=x)
    # factors, standard errors and adjustments
    Vinv <- ginv(t(x) %*% diag(1/(h+d)) %*% x)
    degf <- dim(x)[1] - dim(x)[2]
    vbar <- qfunc(y,d,x)(h) / degf
    Cfct <- sqrt(diag(Vinv))
    Afct <- sqrt(diag(Vinv) *        (vbar))
    Kfct <- sqrt(diag(Vinv) * max(1, (vbar)))
    # critical values
    crt        <- function(qval) {return(qt(qval, df=degf))}
    classic    <- makeConfInts(sgn=sgnf, pst=By, fct=Cfct
                       , crt=qnorm, "unadjusted likelihood")
    knappAdjst <- makeConfInts(sgn=sgnf, pst=By, fct=Afct
                       , crt=crt, "Knapp-Hartung adjustment")
    knappAdhoc <- makeConfInts(sgn=sgnf, pst=By, fct=Kfct
                       , crt=crt, "Knapp-Hartung ad hoc improvement")
    confints   <- rbind(classic, knappAdjst, knappAdhoc)
    confints$h <- names(h) # Don't need this line
    return(confints)
}

#' Interval estimates: For the regression coefficients
#'
#' @param y      study responses.
#' @param d      heteroscedasticity.
#' @param h_dat  data frame of tau estimates.
#' @param x      design matrix.
#' @param sgnf   significance levels.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#' bcg_h <- hEstimates(y=bcg_y, d=bcg_d, x=bcg_x)
#' sgnf_lev <- c(0.01, 0.025, 0.05, 0.01)
#'
#' intervalEstimates(y=bcg_y, d=bcg_d, h_dat=bcg_h, x=bcg_x, sgnf=0.025)
#' intervalEstimates(y=bcg_y, d=bcg_d, h_dat=bcg_h, x=bcg_x,
#'   sgnf=sgnf_lev)
#' @export
intervalEstimates <- function (  y     # study responses
                               , d     # heteroscedasticity
                               , h_dat # data frame of tau estimates
                               , x     # design matrix
                               , sgnf  # significance levels
                               ) {
    ## Confidence intervals based on first order limit theory.
    func <- function (r) {
        if (is.na(r$h)) data.frame() else {
            intervalEstimates_OneSgnf(y,d,r$h,x,sgnf)}
    }
    return(ddply(h_dat, "type", func))
}


#' Inference: Based on methods of moments and
#' maximum likelihood.
#'
#' Calculates the so called Q-profiling confidence interval for
#' the heterogeneity for data following a random effects meta
#' regression model.
#' @param y k-vector of study responses.
#' @param d k-vector of heteroscedasticity.
#' @param x design k-p-matrix.
#' @param sgnf significance levels.
#' @return A data frame containing the bounds of the interval estimate.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#' sgnf_lev <- c(0.01, 0.025, 0.05, 0.01)
#'
#' set.seed(865287113) # for reproducibility
#'
#' hConfidence(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=0.025)
#' hConfidence(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=sgnf_lev)
#' @export
hConfidence <- function (  y    # k-vector of study responses
                         , d    # k-vector of heteroscedasticity
                         , x    # design k-p-matrix
                         , sgnf      # sigificance levels
                         ) {
    bounds <- pfunc(y=y,d=d,x=x)(qchisq(  c(1-sgnf/2,sgnf/2)
                                        , dim(x)[1] - dim(x)[2]))
    return(data.frame(type="q-profiling"
                      , confidence=1-sgnf
                      , lower=bounds[1:length(sgnf)]
                      , upper=bounds[-(1:length(sgnf))]))
}

#' Inference: Based on methods of moments and
#' maximum likelihood.
#'
#' Calculates common statistics for point and confidence interval
#' estimates for the heterogeneity and the regression coefficients
#' of the random effects meta regression model based on the given
#' data.
#'
#' @param y    k-vector of study responses.
#' @param d    k-vector of heteroscedasticity.
#' @param x    design k-p-matrix.
#' @param sgnf significance levels.
#' @return The same return type as the skeleton 'metagenEmpty()'.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#' sgnf_lev <- c(0.01, 0.025, 0.05, 0.01)
#'
#' set.seed(865287113) # for reproducibility
#'
#' c1 <- metareg(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=0.025)
#' c2 <- metareg(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=sgnf_lev)
#'
#' # When performing a meta analysis, provide the function
#' # with a vector of 1s.
#' if (!all(names(c1) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(c2) == names(metagenEmpty()))) stop("Name clash")
#' @export
metareg <- function (  y    # k-vector of study responses
                     , d    # k-vector of heteroscedasticity
                     , x    # design k-p-matrix
                     , sgnf # significance levels
                     ) {
    if (is.vector(x)) {x <- as.matrix(x, ncol=1)}
    h_dat  <- hEstimates(y,d,x)
    # pointr may include an empty data frame for NA h-estimates.
    pointr <- regressionEstimates(y,d,h_dat,x)
    # confr may include an empty data frame for NA h-estimates.
    confr  <- intervalEstimates(y,d,h_dat,x,sgnf)
    confh  <- hConfidence(y=y,d=d,x=x,sgnf=sgnf)
    return(list(pointh=h_dat, confh=confh, pointr=pointr, confr=confr))
}

#' Inference: Empty skeleton
#'
#' Returns an empty skeleton that has
#' the same return type as any other
#' 'metagenSOMETHING' function.
#'
#' @examples
#' metagenEmpty()
#' @export
metagenEmpty <- function () {
    return(list(  pointh=data.frame()
                , confh=data.frame()
                , pointr=data.frame()
                , confr=data.frame())
    )
}

#' Inference: Analysis of the data set
#'
#' Runs all implemented methods and combines them
#' in a neat summary.
#'
#' @param y       k-vector of responses.
#' @param d       k-vector of heteroscedasticities.
#' @param x       design k-p-matrix.
#' @param sgnf    vector of significance levels.
#' @param s       k-vector of study responses. Default is NULL. If
#' 'adjusted=TRUE', this value needs to be given.
#' @param n      draws from the pivotal distribution.
#' @param method Default is 'list("univariate", "multivariate")'.
#' @param adjusted : TRUE or FALSE
#' @return The same return type as the skeleton 'metagenEmpty()'.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#' sgnf_lev <- c(0.01, 0.025, 0.05, 0.01)
#'
#' set.seed(865287113) # for reproducibility
#'
#' # Runs a standard analysis, use n=1000 in an actual
#' # analysis instead!
#' m1 <- metagen(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=0.025, n=50)
#' m2 <- metagen(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=sgnf_lev, n=50)
#'
#' # Runs the methods based on generalised principles via an
#' # adjustment for the unknown heteroscedasticity.  Use
#' # n=1000 in an actual analysis instead!!
#' bcg_s <- bcg$size
#' m3 <- metagen(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=0.025, s=bcg_s, n=50,
#'   adj=TRUE)
#' m4 <- metagen(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=sgnf_lev, s=bcg_s,
#'   n=50, adj=TRUE)
#'
#' if (!all(names(m1) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(m2) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(m3) == names(metagenEmpty()))) stop("Name clash")
#' if (!all(names(m4) == names(metagenEmpty()))) stop("Name clash")
#' @export
metagen <- function (  y      # k-vector of responses.
                     , d      # k-vector of heteroscedasticities.
                     , x      # design k-p-matrix.
                     , sgnf   # vector of significance levels
                     , s=NULL # k-vector of study responses
                     , n      # draws from the pivotal distribution.
                     , method=list("univariate", "multivariate")
                     , adjusted=FALSE # TRUE or FALSE
                     ) {
    cls    <- metareg(y=y, d=d, x=x, sgnf=sgnf)
    gen    <- metagenGeneralised(y=y, d=d, x=x, sgnf=sgnf, n=n,
                                 adjusted=FALSE)
    genAdj <- if (adjusted) {
        metagenGeneralised(y=y, d=d, x=x, sgnf=sgnf, s=s, n=n,
                           adjusted=TRUE)
    } else {metagenEmpty()}

    return(Map(rbind, cls, gen, genAdj))
}
