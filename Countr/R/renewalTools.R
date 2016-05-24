.distDescr <- list(
    ##        "custom" = {
    ##            customPars <- list(...)[[1]]
    ##            parNames <- customPars[["parNames"]]
    ##            linkDefaultNames <- rep("identity", length(parNames))
    ##            names(linkDefaultNames) <- parNames
    ##        },
    "weibull" = list(
        parNames         = c("scale", "shape"),
        linkDefaultNames = c(scale = "log", shape = "log")
    ),

    "gamma" = list(
        parNames         = c("rate", "shape"),
        linkDefaultNames = c(rate = "log", shape = "log")
    ),

    "gengamma" = list(
        parNames         = c("mu", "sigma", "Q"),
        linkDefaultNames = c(mu = "identity", sigma = "log", Q = "identity")
    ),

    "burr" = list(
        parNames         = c("scale", "shape1", "shape2"),
        linkDefaultNames = c(scale = "log", shape1 = "log", shape2 = "log")
    )
)



.getParNames <- function(dist, ...) {
    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               customPars[["parNames"]]
           },
           ## default
           .distDescr[[dist]]$parNames
           )
}

.getDefaultLinkNames <- function(dist, parNames) {
    switch(dist,
           "custom" = {
               linkDefaultNames <- rep("identity", length(parNames))
               names(linkDefaultNames) <- parNames
               linkDefaultNames
           },
           ## default
           .distDescr[[dist]]$linkDefaultNames
           )
}

## inverse links - no need to compute them every time.
.inverseLinks <- list(
    "log"      = function(x) exp(x),
    "cauchit"  = function(x) cauchit(x, inverse = TRUE),
    "cloglog"  = function(x) cloglog(x, inverse = TRUE),
    "probit"   = function(x) probit(x, inverse = TRUE),
    "logit"    = function(x) logit(x, inverse = TRUE),
    "identity" = function(x) x
)
## set attribute "functionName" for the functions in .inverseLinks
##  (note the non-local assignment)
lapply(names(.inverseLinks),
       function(link) attr(.inverseLinks[[link]], "functionName") <<- link
       )

## compute the link function inverse.
.computeInverseLink <- function(link = c("log", "cauchit", "cloglog",
                                    "probit", "logit", "identity")) {
    link <- match.arg(link)
    obj <- switch(link,
                  "log" = {function(x) exp(x)},
                  "cauchit" = {function(x) cauchit(x, inverse = TRUE)},
                  "cloglog" = {function(x) cloglog(x, inverse = TRUE)},
                  "probit" = {function(x) probit(x, inverse = TRUE)},
                  "logit" = {function(x) logit(x, inverse = TRUE)},
                  "identity" = {function(x) x}
                  )
    attr(obj, "functionName") <- link
    obj
}

## TODO: remove after testing and replace the calls to  .computeInverseLink() accordingly.
stopifnot(all(sapply(names(.inverseLinks),
    function(link)
        identical(body(.inverseLinks[[link]]), body(.computeInverseLink(link))) &&
        identical(args(.inverseLinks[[link]]), args(.computeInverseLink(link))) &&
        identical(attr(.inverseLinks[[link]], "functionName"),
                  attr(.computeInverseLink(link), "functionName"))
    )))




## create a string with link function information
.summarizeLinkInformation <- function(linkObj) {
    res <- character()
    for (i in seq(along = linkObj))
        res[i] <- paste0(names(linkObj)[i], ": link ",
                         attr(linkObj[[i]], "functionName"))
    paste0(res, collapse = ", ")
}



#' Creates the renewal control list
#'
#' Creates the renewal control list used by \code{renewal}
#'
#' The function usues user passed inputs, checks them and return an
#' appropriate lists that is used inside \code{renewal} by the optimization
#' rourine \code{optimx} among others
#' @param method character one of the optimization method accepted by
#' \code{optimx}
#' @param maxit numeric the maximum number of iterations in the
#' optimization routine.
#' @param trace Non-negative integer. Should tracing information by
#' printed to the screen.
#' @param start (named) numeric, vector of starting values.
#' @param kkt locical should the Kuhn, Karush, Tucker optimality
#' conditions be tested ? default to \code{FALSE} to avoid numerical
#' hessian computation.
#' @param ... TODO
#' @return a list with the control parameters used by \code{renewal}
#' @export
renewal.control <- function(method = "nlminb", maxit = 1000, trace = 1,
                            start = NULL, kkt = FALSE, ...) {

    rval <- list(method = method, maxit = maxit, trace = trace, kkt = kkt,
                 start = start)
    rval <- c(rval, list(...))
    rval$maximize <- TRUE

    rval
}

#' Creates the convolution inputs setting
#'
#' Check and creates the convolution inputs list
#'
#' @param convPars  a list of convolution parameters arguments with slots
#' \code{nsteps}, \code{extrap} and \code{convMethod}.
#' See \code{dCount_conv_bi}. If NULL, default parameters will be applied.
#' @param dist TODO
#' @return list convolution inputs.
#' @export
renewal.convPars <- function(convPars, dist) {
    extrap <- ifelse(!is.null(convPars$extrap), convPars$extrap,
                     ifelse(dist == "custom", FALSE, TRUE)
                     )
    nsteps <- ifelse(!is.null(convPars$nsteps),
                     convPars$nsteps,
                     ifelse(extrap, 50, 150))
    convMethod <- ifelse(!is.null(convPars$convMethod),
                         convPars$convMethod, "dePril")

    list(nsteps = nsteps, extrap = extrap, convMethod = convMethod)
}

## check custom-parameters
.checkcustomPars <- function(customPars, extrap) {
    ## check that par names exists
    if (is.null(customPars$parNames))
        stop("parNames should be provided in customPars !")
    else {
        if (length(customPars$parNames) < 1)
            stop("parNames must be at least of length 1 !")
        if (!is.character(customPars$parNames))
            stop("parNames must be a character string !")
    }

    ## check survivalFct
    if (is.null(customPars$survivalFct))
        stop("survivalFct should be provided in customPars !")
    else if (class(customPars$survivalFct) != "function")
        stop("survivalFct must be a function object !")

    if (extrap) {
        ## check extrapolFct
        if (is.null(customPars$extrapolFct))
            stop("extrapolFct should be provided in customPars !")
        else if (class(customPars$extrapolFct) != "function")
            stop("extrapolFct must be a function object !")
    }

    customPars
}


#' Creates the series expansion inputs setting
#'
#' Check and creates the series expansion inputs list
#'
#' @param seriesPars list series expansion input parameters with slots
#' \code{terms} (number of terms in the series expansion),
#' \code{iter} (number of iteration in the accelerated series expansion
#' algorithm) and \code{eps} (tolerance in the accelerated series expansion
#' algorithm), Only used if \code{dist = "weibull"} and
#' \code{weiMethod = c("series_mat", "series_acc")}.
#' @param long TODO
#' @return list series expansion inputs.
#' @export
renewal.seriesPars <- function(seriesPars, long = FALSE) {
    terms <- ifelse(!is.null(seriesPars$terms), seriesPars$terms,
                    ifelse(long, 100, 50))
    iter <- ifelse(!is.null(seriesPars$iter), seriesPars$iter, 300)
    eps <- ifelse(!is.null(seriesPars$eps), seriesPars$eps, 1e-10)

    list(terms = terms, iter = iter, eps = eps)
}

#' Check weibull computation algorithm
#'
#' Check weibull computation algorithm
#'
#' @param weiMethod character weibull method desired.
#' @return a valid weibull computation method.
#' @export
renewal.weiMethod <- function(weiMethod) {
    if (is.null(weiMethod))
        weiMethod <- "series_acc"
    else if (!weiMethod %in% c("series_acc", "series_mat",
                               "conv_direct", "conv_naive", "conv_dePril")) {
        weiMethod <- "series_acc"
        warning(paste(weiMethod, "is not an accepted method for weibull dist!",
                      "accelerated series will be used !"))
    }

    weiMethod
}


## compute the link function inverse.
#' @importFrom VGAM cauchit
#' @importFrom VGAM cloglog
#' @importFrom VGAM probit
#' @importFrom VGAM logit
.computeInverseLink <- function(link = c("log", "cauchit", "cloglog",
                                    "probit", "logit", "identity")) {
    link <- match.arg(link)
    obj <- switch(link,
                  "log" = {function(x) exp(x)},
                  "cauchit" = {function(x) cauchit(x, inverse = TRUE)},
                  "cloglog" = {function(x) cloglog(x, inverse = TRUE)},
                  "probit" = {function(x) probit(x, inverse = TRUE)},
                  "logit" = {function(x) logit(x, inverse = TRUE)},
                  "identity" = {function(x) x}
                  )
    attr(obj, "functionName") <- link
    obj
}

## create a string with link function information
.summarizeLinkInformation <- function(linkObj) {
    res <- character()
    for (i in seq(along = linkObj))
        res[i] <- paste0(names(linkObj)[i], ": link ",
                         attr(linkObj[[i]], "functionName"))
    paste0(res, collapse = ", ")
}

## apply (inverse) link function to the standard errors of the parameters
.transformSE <- function(se, linkList) {
    nmPars <- names(linkList)
    ## i=1 always transform
    par1 <- nmPars[1]
    key1 <- paste0(par1, "_")
    ll1 <- which(grepl(key1, names(se)))
    se[ll1] <- linkList[[1]](se[ll1])

    for (i in 2:length(nmPars)) {
        pari <- nmPars[i]
        keyi <- paste0(pari, "_")
        lli <- which(grepl(keyi, names(se)))
        if (length(lli) > 1)
            se[lli] <- linkList[[i]](se[lli])
        else if (length(lli) == 1) {
            if (nchar(names(se)[lli]) > nchar(keyi)) ##
                se[lli] <- linkList[[i]](se[lli])
        }
    }
    se
}


.getDistParNames <- function(dist = c("weibull", "weibullgam",
                                 "gamma", "gengamma", "burr")) {
    dist <- match.arg(dist)

    switch(dist,
           "weibull" = {
               parNames <- c("scale", "shape")
           },
           "gamma" = {
               parNames <- c("rate", "shape")
           },
           "gengamma" = {
               parNames <- c("mu", "sigma", "Q")
           },
           "burr" = {
               parNames <- c("scale", "shape1", "shape2")
           },
           "weibullgam" = {
               parNames <- c("scale", "shape", "shapeGam", "scaleGam")
           }
           )

    parNames
}

.getLinkDefaultNames <- function(dist = c("weibull", "weibullgam",
                                     "gamma", "gengamma", "burr")) {
    dist <- match.arg(dist)
    switch(dist,
           "weibull" = {
               linkDefaultNames <- c(scale = "log", shape = "log")
           },
           "gamma" = {
               linkDefaultNames <- c(rate = "log", shape = "log")
           },
           "gengamma" = {
               linkDefaultNames <- c(mu = "identity", sigma = "log",
                                     Q = "identity")
           },
           "burr" = {
               linkDefaultNames <- c(scale = "log", shape1 = "log",
                                     shape2 = "log")
           },
           "weibullgam" = {
               linkDefaultNames <- c(scale = "log",
                                     shape = "identity",
                                     shapeGam = "identity",
                                     scaleGam = "identity")
           }
           )
    linkDefaultNames
}

## get the inverse link function for each distribution
.getLinkList <- function(dist, link, ...) {

    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
               linkDefaultNames <- rep("identity", length(parNames))
               names(linkDefaultNames) <- parNames
           }, {
               parNames <- .getDistParNames(dist)
               linkDefaultNames <- .getLinkDefaultNames(dist)
           }
           )

    linkRes <- list()
    for (i in seq(along = parNames)) {
        pari <- parNames[i]
        linki <- link[[pari]]
        if(is.null(linki))
            linki <- .computeInverseLink(linkDefaultNames[pari])
        else
            linki <- .computeInverseLink(linki)

        linkRes[[pari]] <- linki
    }

    class(linkRes) <- "InverseLink"
    linkRes
}

## model matrix returned as a list
.getModelMatrix <- function(formula, dist, mf, anc, ...) {
    X <- model.matrix(formula, data = mf, rhs = 1)
    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
           },
           {
               parNames <- .getDistParNames(dist)
           }
           )

    modelMatRes <- list()
    modelMatRes[[parNames[1]]] <- X
    if (length(parNames) > 1) {
        for (i in 2:length(parNames)) {
            pari <- parNames[i]
            if (pari %in% names(anc)) {
                fi <- Formula(anc[[pari]])
                Xi <- model.matrix(fi, data = mf, rhs = 1)
                modelMatRes[[pari]] <- Xi
            }
        }
    }

    modelMatRes
}

## check initial values
.checkInitialValues <- function(dist, start, modelMatrixList, weights, Y,
                                anc, ...) {
    switch(dist,
           "weibullgam" ={
               X <- modelMatrixList[["scale"]]
               nm <- c(paste0("scale_", colnames(X)),
                       "shape_", "shapeGam_", "scaleGam_")

               valid <- FALSE
               if (!is.null(start)) {
                   if(all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unamed initial values found !",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   IV <- glm.fit(X, Y, family = poisson(), weights = weights)
                   ## change parametrization in terms of 1/r and 1/alpha (the gamma pars)
                   start <- c(IV$coefficients, 1, 0.061, 0.061 * 3)
                   names(start) <- nm
               }

           },
           "weibull" = {
               X <- modelMatrixList[["scale"]]
               nm <- paste0("scale_", colnames(X))
               if ("shape" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("shape_",
                                      colnames(modelMatrixList[["shape"]])
                                      )
                           )
               else
                   nm <- c(nm, "shape_")

               valid <- FALSE
               if (!is.null(start)) {
                   if(all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unamed initial values found !",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (!"shape" %in% names(modelMatrixList)) {
                       IV <- glm.fit(X, Y, family = poisson(), weights = weights)
                       start <- c(IV$coefficients, 1)
                       names(start) <- nm
                   } else
                       stop(paste("initial values should be provided in control",
                                  "when regression on ancillary parameters !"))
               }
           },
           "gamma" = {
               X <- modelMatrixList[["rate"]]
               nm <- paste0("rate_", colnames(X))
               if ("shape" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("shape_",
                                      colnames(modelMatrixList[["shape"]])
                                      )
                           )
               else
                   nm <- c(nm, "shape_")

               valid <- FALSE
               if (!is.null(start)) {
                   if(all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unamed initial values found !",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (!"shape" %in% names(modelMatrixList)) {
                       IV <- glm.fit(X, Y, family = poisson(), weights = weights)
                       start <- c(IV$coefficients, 1)
                       names(start) <- nm
                   } else
                       stop(paste("initial values should be provided in control",
                                  "when regression on ancillary parameters !"))
               }
           },
           "gengamma" = {
               X <- modelMatrixList[["mu"]]
               nm <- paste0("mu_", colnames(X))

               if ("sigma" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("sigma_",
                                      colnames(modelMatrixList[["sigma"]])))
               else
                   nm <- c(nm, "sigma_")

               if ("Q" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("Q_",
                                      colnames(modelMatrixList[["Q"]])))
               else
                   nm <- c(nm, "Q_")

               valid <- FALSE
               if (!is.null(start)) {
                   if (all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unamed initial values found !",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (!any(c("sigma" %in%  names(modelMatrixList),
                              "Q" %in%  names(modelMatrixList))
                            )
                       ) {
                       IV <- glm.fit(X, Y, family = poisson(), weights = weights)
                       start <- c(IV$coefficients, 1, 1)
                       names(start) <- nm
                   } else
                       stop(paste("initial values should be provided in control",
                                  "when regression on ancillary parameters !"))
               }
           },
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]

               location <- parNames[1]
               X <- modelMatrixList[[location]]
               nm <- paste0(location, "_", colnames(X))

                if (length(parNames) > 1) {
                   for (i in 2:length(parNames)) {
                       pari <- parNames[i]
                       if (pari %in% names(anc)) {
                           nm <- c(nm, paste0(pari, "_",
                                              colnames(modelMatrixList[[pari]])
                                              )
                                   )
                       } else
                           nm <- c(nm, paste0(pari, "_"))
                   }
               }

               valid <- FALSE
               if (!is.null(start)) {
                   if (all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unamed initial values found !",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid)
                   stop(paste("initial values should be provided in control",
                              "when custom survival functions are passed !"))
           },
           "burr" = {
               X <- modelMatrixList[["scale"]]
               nm <- paste0("scale_", colnames(X))

               if ("shape1" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("shape1_",
                                      colnames(modelMatrixList[["shape1"]])))
               else
                   nm <- c(nm, "shape1_")

               if ("shape2" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("shape2_",
                                      colnames(modelMatrixList[["shape2"]])))
               else
                   nm <- c(nm, "shape2_")

               valid <- FALSE
               if (!is.null(start)) {
                   if (all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unamed initial values found !",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (!any(c("shape1" %in%  names(modelMatrixList),
                              "shape2" %in%  names(modelMatrixList))
                            )
                       ) {
                       IV <- glm.fit(X, Y, family = poisson(), weights = weights)
                       start <- c(IV$coefficients, 1, 3)
                       names(start) <- nm
                   } else
                       stop(paste("initial values should be provided in control",
                                  "when regression on ancillary parameters !"))
               }
           }
           )
    start
}

## return the loglikelihood (if Ev = FALSE) or the expected value of the
## count if Ev = TRUE
.objectiveFunction <- function(params, dist, modelMatrixList, linkList,
                               time, convPars, Y = NULL, weights = NULL,
                               Ev = FALSE, summa = TRUE,
                               seriesPars = NULL, weiMethod = NULL,
                               ...) {

    ## check if Y is provided when Ev = TRUE
    if (!Ev) {
        if(is.null(Y))
            stop("Y should be provided when Ev is FALSE!")
        else {
            ## weights
            if (is.null(weights))
                weights <- rep(1, length(Y))
        }
    }

    ## get xMax,
    dotPars <- list(...)
    if (!is.null(dotPars$xMax))
        xMax <- dotPars$xMax
    else if (!is.null(Y))
        xMax <- max(max(Y), 15)
    else {
        xMax <- 15
        print("xMax not provided! 15 will be used by default !")
    }

    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
               distParValues <- list()
               locationPar <- parNames[1]
               linkLocation <- linkList[[locationPar]]
               X <- modelMatrixList[[locationPar]]

               ## we assume regression for the first parameter (location)
               if(is.null(names(params)))
                   location <-
                       as.vector(linkLocation(X %*% params[1:ncol(X)]))
               else
                   location <-
                       as.vector(linkLocation(X %*%
                                              params[grepl(paste0(locationPar,
                                                                  "_"),
                                                           names(params))]
                                              )
                                 )

               distParValues[[locationPar]] <- location
               ix <- ncol(X)
               for (i in 2:length(parNames)) {
                   pari <- parNames[i]
                   if (pari %in% names(modelMatrixList)) {
                       linki <- linkList[[pari]]
                       Xi <- modelMatrixList[[pari]]
                       if(is.null(names(params)))
                           distParValues[[pari]] <-
                               as.vector(
                                   linki(Xi %*%
                                         params[(ix + 1):
                                                (ix + ncol(Xi))
                                                ]
                                         )
                                   )
                       else
                           distParValues[[pari]] <-
                               as.vector(
                                   linki(Xi %*%
                                         params[grepl(paste0(pari, "_"),
                                                      names(params))]
                                         )
                                   )

                       ix <- ix + ncol(Xi)
                   } else {
                       distParValues[[pari]] <-
                           as.numeric(rep(params[ix + 1], nrow(X)))
                       ix <- ix + 1
                   }
               }

               ## seq along Y
               seky <- seq(along = Y)
               ## extraction tools
               survFct <- customPars[["survivalFct"]]
               ## create distPars in the appropriate format
               .fct0 <- function(xlist, i = 1) xlist[[i]]
               .fcti <- function(i) lapply(distParValues, .fct0, i = i)
               distPars <- lapply(seky, .fcti)

               if (convPars$extrap) {
                   .getextrapolPars <- customPars[["extrapolFct"]]
                   .fctExtrap <- function(i) .getextrapolPars(distPars[[i]])
               } else
                   .fctExtrap <- function(i) c(2, 2) ## not used

               extrapPrasList <- lapply(seky, .fctExtrap)

               if (Ev) { ## expected value and variance
                   .evfct <- function(i)
                       evCount_conv_user(xmax = xMax, survR = survFct,
                                         distPars = distPars[[i]],
                                         extrapolPars = .fctExtrap(i),
                                         method = convPars$convMethod,
                                         nsteps = convPars$nsteps,
                                         time = time,
                                         extrap = convPars$extrap)

                   return(lapply(seky, .evfct))

               } else { ## loglikelihood or probability value
                   if (summa)
                       return(dCount_conv_loglik_user(x = Y,
                                                      survR = survFct,
                                                      distPars = distPars,
                                                      extrapolPars =
                                                      extrapPrasList,
                                                      nsteps = convPars$nsteps,
                                                      extrap = convPars$extrap,
                                                      method =
                                                      convPars$convMethod,
                                                      time = time,
                                                      weights = weights)
                              )
                   else {
                       .vectoGetProbs <- function(i) {
                           dCount_conv_user(x = Y[i], survR = survFct,
                                            distPars = distPars[[i]],
                                            extrapolPars = .fctExtrap(i),
                                            nsteps = convPars$nsteps,
                                            extrap = convPars$extrap,
                                            method = convPars$convMethod,
                                            time = time)
                           }

                       return(as.numeric(sapply(seky, .vectoGetProbs)))
                   }
               }
           },
           "weibull" = { ## we do weibull separately on purpose
               if (is.null(seriesPars))
                   seriesPars <- renewal.seriesPars(seriesPars)
               if (is.null(weiMethod))
                   weiMethod <- renewal.weiMethod(weiMethod)
               ## ------------- scale parameter -----------------------------
               linkScale <- linkList[["scale"]]
               X <- modelMatrixList[["scale"]]
               ## check parameters name
               if(is.null(names(params)))
                   scale <-
                       as.vector(linkScale(X %*% params[1:ncol(X)]))
               else
                   scale <-
                       as.vector(linkScale(X %*%
                                           params[grepl("scale_",
                                                        names(params))]
                                           )
                                 )

               ## ------------- shape parameter -----------------------------
               if ("shape" %in% names(modelMatrixList)) {
                   linkShape <- linkList[["shape"]]
                   Xshape <- modelMatrixList[["shape"]]
                   if(is.null(names(params)))
                       shape <- as.vector(
                           linkShape(Xshape %*%
                                     params[(ncol(X) + 1):length(params)]
                                     )
                           )
                   else
                       shape <- as.vector(
                           linkShape(Xshape %*%
                                     params[grepl("shape_", names(params))]
                                     )
                           )
               } else
                   shape <- rep(params[length(params)], nrow(X))

               if (Ev) {
                   .fct <- function(i)
                       evWeibullCount(xmax = xMax, shape = shape[i],
                                      scale = scale[i],
                                      method = weiMethod, time = time,
                                      conv_steps = convPars$nsteps,
                                      conv_extrap = convPars$extrap,
                                      series_terms = seriesPars$terms,
                                      series_acc_niter = seriesPars$iter,
                                      series_acc_eps = seriesPars$eps
                                      )

                   return(lapply(1:nrow(X), .fct))

               } else {
                   if (summa) {
                       return(
                           dWeibullCount_loglik(
                               x = Y, shape = shape, scale = scale,
                               conv_steps = convPars$nsteps,
                               conv_extrap = convPars$extrap,
                               method = weiMethod,
                               time = time, weights = weights,
                               series_terms = seriesPars$terms,
                               series_acc_niter = seriesPars$iter,
                               series_acc_eps = seriesPars$eps)
                              )
                   } else
                       .vectoGetProbs <- function(i) {
                           dWeibullCount(
                               x = Y[i], shape = shape[i], scale = scale[i],
                               conv_steps = convPars$nsteps,
                               conv_extrap = convPars$extrap,
                               method = weiMethod, time = time,
                               series_terms = seriesPars$terms,
                               series_acc_niter = seriesPars$iter,
                               series_acc_eps = seriesPars$eps)
                       }

                   return(as.numeric(sapply(seq(along = Y),
                                            .vectoGetProbs))
                          )
               }
           },
           "weibullgam" = { ## we do weibull-gamma separately as well
               ## ------------- scale parameter -----------------------------
               X <- modelMatrixList[["scale"]]
               ## check parameters name
               if(is.null(names(params))) {
                   scale <- as.vector(params[1:ncol(X)])
                   shape <- params[ncol(X) + 1]
                   shapeGam <- params[ncol(X) + 2]
                   scaleGam <- params[ncol(X) + 3]
               } else {
                   scale <- params[grepl("scale_", names(params))]
                   shape <- params["shape_"]
                   shapeGam <- params["shapeGam_"]
                   scaleGam <- params["scaleGam_"]
               }

               ## Note that we change parametrization here in terms of
               ## 1/r and 1/alpha

               if (Ev) {
                   .fct <- function(i) {
                       evWeibullgammaCount(xmax = xMax, shape = shape,
                                           shapeGam = 1.0 / shapeGam,
                                           scaleGam = 1.0 / scaleGam,
                                           Xcovar = X[i, ,drop = FALSE],
                                           beta = scale,
                                           method = weiMethod, time = time,
                                           series_terms = seriesPars$terms,
                                           series_acc_niter = seriesPars$iter,
                                           series_acc_eps = seriesPars$eps
                                      )
                   }
                   return(lapply(1:nrow(X), .fct))

               } else {
                   if (summa) {
                       return(
                           dWeibullgammaCount_loglik(
                               x = Y, shape = shape,
                               shapeGam = 1.0 / shapeGam, scaleGam = 1.0 / scaleGam,
                               Xcovar = X, beta = scale,
                               method = weiMethod,
                               time = time, weights = weights,
                               series_terms = seriesPars$terms,
                               series_acc_niter = seriesPars$iter,
                               series_acc_eps = seriesPars$eps)
                              )
                   } else
                       return(
                           dWeibullgammaCount(
                               x = Y, shape = shape,
                               shapeGam = 1.0 / shapeGam,
                               scaleGam = 1.0 / scaleGam,
                               Xcovar = X, beta = scale,
                               method = weiMethod, time = time,
                               series_terms = seriesPars$terms,
                               series_acc_niter = seriesPars$iter,
                               series_acc_eps = seriesPars$eps)
                              )
               }
           },
           { ## default
               parNames <- .getDistParNames(dist)
               distParValues <- list()
               locationPar <- parNames[1]
               linkLocation <- linkList[[locationPar]]
               X <- modelMatrixList[[locationPar]]

               ## we assume regression for the first parameter (location)
               if(is.null(names(params)))
                   location <-
                       as.vector(linkLocation(X %*% params[1:ncol(X)]))
               else
                   location <-
                       as.vector(linkLocation(X %*%
                                              params[grepl(paste0(locationPar,
                                                                  "_"),
                                                           names(params))]
                                              )
                                 )

               distParValues[[locationPar]] <- location
               ix <- ncol(X)
               for (i in 2:length(parNames)) {
                   pari <- parNames[i]
                   if (pari %in% names(modelMatrixList)) {
                       linki <- linkList[[pari]]
                       Xi <- modelMatrixList[[pari]]
                       if(is.null(names(params)))
                           distParValues[[pari]] <-
                               as.vector(
                                   linki(Xi %*%
                                         params[(ix + 1):
                                                (ix + ncol(Xi))
                                                ]
                                         )
                                   )
                       else
                           distParValues[[pari]] <-
                               as.vector(
                                   linki(Xi %*%
                                         params[grepl(paste0(pari, "_"),
                                                      names(params))]
                                         )
                                   )

                       ix <- ix + ncol(Xi)
                   } else {
                       distParValues[[pari]] <-
                           as.numeric(rep(params[ix + 1], nrow(X)))
                       ix <- ix + 1
                   }
               }

               ## seq along Y
               seky <- seq(along = Y)
               ## create distPars in the appropriate format

               .fct0 <- function(xlist, i = 1) xlist[[i]]
               .fcti <- function(i) lapply(distParValues, .fct0, i = i)
               distPars <- lapply(seky, .fcti)

               if (Ev) { ## expected value and variance
                   .evfct <- function(i)
                       evCount_conv_bi(xmax = xMax,
                                       distPars = distPars[[i]],
                                       dist = dist,
                                       method = convPars$convMethod,
                                       nsteps = convPars$nsteps,
                                       time = time,
                                       extrap = convPars$extrap)

                   return(lapply(seky, .evfct))

               } else { ## loglikelihood or probability value
                   if (summa)
                       return(dCount_conv_loglik_bi(x = Y,
                                                    distPars = distPars,
                                                    dist = dist,
                                                    nsteps = convPars$nsteps,
                                                    extrap = convPars$extrap,
                                                    method =
                                                    convPars$convMethod,
                                                    time = time,
                                                    weights = weights)
                              )
                   else {
                       .vectoGetProbs <- function(i) {
                           dCount_conv_bi(x = Y[i],
                                          distPars = distPars[[i]],
                                          dist = dist,
                                          nsteps = convPars$nsteps,
                                          extrap = convPars$extrap,
                                          method = convPars$convMethod,
                                          time = time)
                           }

                       return(as.numeric(sapply(seky, .vectoGetProbs,)))
                   }
               }
           })
}

## extract apprpriate modelData that can be used later to construc se
.modelData <- function(modelMatrixList, dist, ...) {
    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
           }, {
               parNames <- .getDistParNames(dist)
           }
           )

    locationPar <- parNames[1]
    X <- modelMatrixList[[locationPar]]
    nm <- paste0(locationPar, "_", colnames(X))
    colnames(X) <- gsub('\\(Intercept\\)', "", nm)

    if (length(parNames) > 1) {
        for (i in 2:length(parNames)) {
            pari <- parNames[i]
            if (pari %in% names(modelMatrixList)) {
                Xi <- modelMatrixList[[pari]]
                nm <- paste0(pari, "_", colnames(Xi))
                colnames(Xi) <- gsub('\\(Intercept\\)', "", nm)
            } else {
                Xi <- matrix(1, ncol = 1, nrow = nrow(X))
                colnames(Xi) <- paste0(pari, "_")
            }
            X <- cbind(X, Xi)
        }
    }
    X
}

.getPredictionStd <- function(modelMatrixList, vcov, dist, link, ...) {
    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
           }, {
               parNames <- .getDistParNames(dist)
           }
           )

    seList <- list()
    ## do location parameter separately to avoid a call to if
    locationPar <- parNames[1]
    X <- modelMatrixList[[locationPar]]
    nm <- paste0(locationPar, "_", colnames(X))
    colnames(X) <- gsub('\\(Intercept\\)', "", nm)

    ## create standard error for the location parameter
    indX <- which(grepl(locationPar, colnames(X)))
    Xloc <- X[, indX]
    indCov <- which(grepl(locationPar, colnames(X)))
    covloc <- vcov[indCov, indCov]
    seloc <- sqrt(diag(Xloc %*% covloc %*% t(Xloc)))
    seList[[locationPar]] <- link[[locationPar]](seloc)

    ## anc parameters
    if (length(parNames) > 1) {
        for (i in 2:length(parNames)) {
            pari <- parNames[i]
            if (pari %in% names(modelMatrixList)) {
                Xi <- modelMatrixList[[pari]]
                nm <- paste0(pari, "_", colnames(Xi))
                colnames(Xi) <- gsub('\\(Intercept\\)', "", nm)
            } else {
                Xi <- matrix(1, ncol = 1, nrow = nrow(X))
                colnames(Xi) <- paste0(pari, "_")
            }
            indXi <- which(grepl(pari, colnames(Xi)))
            Xi <- Xi[, indXi, drop = FALSE]
            indCov <- which(grepl(pari, colnames(Xi)))
            covi <- vcov[indCov, indCov, drop = FALSE]
            sei <- sqrt(diag(Xi %*% covi %*% t(Xi)))
            seList[[pari]] <- link[[pari]](sei)
        }
    }
    seList
}


## temp function
.extractElem <- function(xList, ind = "ExpectedValue")
    as.numeric(xList[[ind]])

.checkHess <- function(hess, nPars) {
    ans <- TRUE

    if (is.null(hess))
        return(FALSE)

    if (any(is.na(hess)))
        return(FALSE)

    if (!is.matrix(hess))
        return(FALSE)

    if (ncol(hess) != nPars)
        return(FALSE)

    if (nrow(hess) != nPars)
        return(FALSE)

    return(ans)
}

## replace NA/NAN that appears in log-likelihood by logMin due to very small
## probability.
#' @keywords internal
.logNaReplace <- function() {
    log(.Machine$double.xmin)
}
