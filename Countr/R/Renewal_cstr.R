#' Fit renewal count processes regression models
#'
#' Fit renewal regression models for count data via maximum likelihood.
#'
#' \code{renewal} re-uses design and functionality of the basic R tools
#' for fitting regression model (\code{lm}, \code{glm}) and is highly
#' inspired from \code{hurdle()} and \code{zeroinfl()} from package
#' \code{pscl}. Package \code{Formula} was used to handle formulae.
#'
#' If User wants to pass his own inter-arrival
#' distribution i.e, \code{dist = "custom"}, he needs to provides his own
#' initial values as well as the \code{customPars} list as follows:
#' \describe{
#' \item{parNames}{character vector the distribution parameters name.
#' Location parameter should be passed first.}
#' \item{survivalFct}{function object containing the survival function. It
#' should have the signature \code{function(t, distPars)} where $t$ is the
#' point where the survival function is evaluated and distPars is the list
#' of the distribution parameters. It should return a double value.}
#' \item{extrapolFct}{function object computing the extrapolation values
#' (numeric of length 2) from the value of the distribution parameters
#' (in \code{distPars}). It should have the signature
#' \code{function(distPars)} and returns a numeric vector of length 2. Only
#' required if the extrapolation is set to \code{TRUE} in \code{convPars}.}
#' }
#' Some checks will be used to validate the \code{customPars} input but
#' it is user responsability to make sure the different functions have the
#' appropriate signatures.
#' \strong{Note:} The weibull gamma is an experimental version and user should
#' use it with care! It is very sensitive to initials values and there is no
#' guarantee of convergence. It has also been reparametrized in terms of
#' (1/r, 1/alpha, c) instead of (r, alpha, c) where r and alpha are the
#' shape and scale of the gamma distribution and c is the shape of the
#' weibull distribution.
#' @param formula single response formula object.
#' @param data,subset,na.action, arguments controlling formula processing
#' via \code{model.frame}. 
#' @param weights optional numeric vector of weights.
#' @param offset optional numeric vector with an a priori known component
#' to be included in the linear predictor of the count model. Currently
#' not used.
#' @param dist character built-in distribution to be used for the inter-
#' arrival time distribution. Currently built in distribution are \code{weibull},
#' \code{weibullgam}, \code{gamma}, \code{gengamma} (generalized-gamma)
#' and \code{burr}. User can provide his own distribution
#' by setting \code{dist} to \code{custom}.
#' @param anc list (named) of formulae to model regression on ancillary
#' parameters. If \code{NULL}, no regression is modelled. Otherwise, the
#' formulae associated with the (exact) parameter name is used.
#' @param link list (named) of character specifiying the name of the link function
#' to be used in the regression. If \code{NULL}, the canonical
#' link function will be used, i.e, \code{log} if the parameter is supposed
#' to be positive, identity otherwise.
#' @param time numeric time at which the count is observed; default to unit (1)
#' @param convPars a list of convolution parameters argumentswith slots
#' \code{nsteps}, \code{extrap} and \code{convMethod}.
#' See \code{dCount_conv_bi}. If NULL, default parameters will be applied.
#' @param control a list of control arguments specified via
#' \code{renewal.control}.
#' @param customPars list user inputs if \code{dist = "custom"}. See details
#' @param seriesPars list series expansion input parameters with slots
#' \code{terms} (number of terms in the series expansion),
#' \code{iter} (number of iteration in the accelerated series expansion
#' algorithm) and \code{eps} (tolerance in the accelerated series expansion
#' algorithm), Only used if \code{dist = "weibull"} and
#' \code{weiMethod = c("series_mat", "series_acc")}.
#' @param weiMethod character computation method to be used if
#' \code{dist = c("weibull", "weibullgam")}. See \code{dWeibullCount} and
#' \code{dWeibullgammaCount}.
#' @param computeHessian logical should the hessian (and hence the covariance
#' matrix) be computed numerically at the fitted values.
#' @param model,y,x logicals. If \code{TRUE} the  corresponding  components
#' of  the  fit  (model  frame,  response, model matrix) are returned.
#' @param ... arguments passed to \code{renewal.control} in the default setup.
#' @return An \code{S3} object of class "renewal", i.e., a list with components
#' including:
#' \describe{
#' \item{coefficients}{value of the fitted coefficients}
#' \item{residuals}{vector of weighted residuals \eqn{\omega * (observed - fitted)}}
#' \item{fitted.values}{vector of fitted means}
#' \item{optim}{data.frame output of \code{optimx}}
#' \item{method}{optimisation algorithm}
#' \item{control}{the control arguments passed to \code{optimx}}
#' \item{start}{starting values  passed to \code{optimx}}
#' \item{weights}{weights applied if any}
#' \item{n}{number of observation (with weights > 0)}
#' \item{iterations}{number of iterations in the optimisation algorithm}
#' \item{execTime}{duration of the optimisation}
#' \item{loglik}{log-likelihood of the fitted model}
#' \item{df.residual}{residuals degrees of freedom for the fitted model}
#' \item{vcoc}{convariance matrix of all coefficients computed numerically from the
#' hessian at the fitted coefficients (if \code{computeHessian} ois \code{TRUE}).}
#' \item{dist}{name of inter-arrival distribution.}
#' \item{link}{list inverse link function corresponding to each parameter in the
#' inter-arrival distribution}
#' \item{converged}{logical did the optimisation algorithm converged ?}
#' \item{data}{data used to fit the model}
#' \item{formula}{the original formula}
#' \item{call}{the original function call}
#' \item{anc}{list (named) of formulae to model regression on ancillary
#' parameters.}
#' \item{convPars}{convolution inputs used}
#' \item{customPars}{user passed distribution inputs. See details}
#' \item{time}{observed window used. default to 1.0 (see inputs)}
#' \item{model}{the full model frame (if \code{model = TRUE}}
#' \item{y}{the response count vector (if \code{y = TRUE}}
#' \item{x}{the model matrix if \code{x = TRUE}}
#' }
#' @export
#' @importFrom numDeriv hessian
#' @importFrom MASS ginv
#' @import optimx
renewal <- function(formula, data, subset, na.action, weights, offset,
                    dist = c("custom", "weibull", "weibullgam",
                        "gamma", "gengamma", "burr"),  
                    anc = NULL, convPars = NULL, link = NULL, time = 1.0,
                    control = renewal.control(...), customPars = NULL,
                    seriesPars = NULL, weiMethod = NULL, computeHessian = TRUE,
                    model = TRUE, y = TRUE, x = FALSE, ...) {

    dist <- match.arg(dist)
    ## check convolution parameters
    convPars <- renewal.convPars(convPars, dist)
    
    if (dist == "custom")
        customPars <- .checkcustomPars(customPars, convPars$extrap)
    else if (dist == "weibull") {
        seriesPars <- renewal.seriesPars(seriesPars)
        weiMethod <- renewal.weiMethod(weiMethod)
    } else if (dist == "weibullgam") {
        warning(
            "weibullgam should be used with care! no guarantee of convergence !")
        anc <- NULL ## no regression allowed on aux pars
        seriesPars <- renewal.seriesPars(seriesPars, TRUE)
        weigamMethod <- weiMethod
        weigamMethod <- ifelse(is.null(weigamMethod), "series_acc",
                               weigamMethod)
        if (! weigamMethod %in% c("series_acc", "series_mat")) {
            warning(paste(weiMethod,
                          "is not an accepted method for weibullgam dist!",
                          "accelerated series will be used !"))
            weigamMethod <- "series_acc"
        }
        weiMethod <- weigamMethod
    }
    
    ## prepare the formula setting
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    
    f <- Formula(formula)
    mf[[1]] <- as.name("model.frame")
    mf$formula <- f
    mf <- eval(mf, parent.frame())
    ## Not sure about this: copied from hurdle: CHECK
    cl <- match.call()
    
    Y <- model.response(mf)
    n <- length(Y)
    ## stop if a formula with multiple response if passed
    if (is.null(Y))
        stop("muti-response formula not accepted !")
    ## convert negative reponse to zeros
    if (length(Y) < 1) 
        stop("empty model")
    if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001))))) 
        warning(paste("invalid dependent variable,",
                      "non-integer values!",
                      "will be transformed")
                )
    Y <- as.integer(round(Y + 0.001))
    if (any(Y < 0)) 
        stop("invalid dependent variable, negative counts")

    ## extract weights and reshape them
    weights <- model.weights(mf)
    if (is.null(weights)) 
        weights <- 1
    if (length(weights) == 1) 
        weights <- rep.int(weights, n)
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)

    ## get the inverse link function
    if (class(link) == "InverseLink")
        linkList <- link
    else
        linkList <- .getLinkList(dist, link, customPars)
    

    ## get model matrices
    modelMatrixList <- .getModelMatrix(f, dist, mf, anc, customPars)
    
    ## create objective function
    ## ============== countDist to pass to optim ===============================
    countDist <- function(params) {
        .objectiveFunction(params, dist, modelMatrixList,
                           linkList, time, convPars, Y, weights,
                           Ev = FALSE, summa = TRUE, seriesPars, weiMethod,
                           customPars)
    }
    
    ## check initilas values
    start <- control$start
    start <- .checkInitialValues(dist, start, modelMatrixList, weights, Y,
                                 anc, customPars)
    nmPars <- gsub('\\(Intercept\\)', "", names(start))
    ## run optimazation routine
    
    method <- control$method
    hessian <- ifelse(is.null(control$hessian), FALSE, control$hessian)
    control$method <- control$start <- control$hessian <- NULL

    if (control$trace)
        print("calling optimx() for parameter estimation by ML ...")
 
    fitCount <- optimx(par = start, fn = countDist, method = method,
                       hessian = hessian, control = control)

    fitCount <- fitCount[1, ]
    ## coefficients
    coefs <- as.numeric(coef(fitCount))
    names(coefs) <- nmPars

    
    ## variance-covariance matrix
    if (computeHessian) {
        hess <- try(attr(fitCount, "details")[method, "nhatend"][[1]])
        if (!.checkHess(hess, length(coefs))) {
            if (control$trace)
                print("computing a numerical approximation to the Hessian ...")
            hess <- numDeriv::hessian(countDist, coefs)
        }
    
        varCovarcount <- try(-solve(hess))
        if(inherits(varCovarcount, "try-error")) {
            varCovarcount <- Matrix::nearPD(-ginv(hess))
            warning(paste("variance-covariance matrix was computed",
                          "by smoothing the genralized inverse hessian !"))
        }
        
        dimnames(varCovarcount) <- list(nmPars, nmPars)
    } else
        varCovarcount <- matrix()
    
    ## residuals (Pearson)
    resTemp <-  .objectiveFunction(coefs, dist, modelMatrixList,
                                   linkList, time, convPars, Y, weights,
                                   TRUE, FALSE, seriesPars, weiMethod,
                                   customPars)
        
    Yhat <- sapply(resTemp, .extractElem, ind = "ExpectedValue")
    wi <- sapply(resTemp, .extractElem, ind = "Variance")
    res <- sqrt(weights) * (Y - Yhat)
    ## number of observations
    nobs <- sum(weights > 0)

    ## value to be returned
    rval <- list(
        coefficients = coefs, residuals = res, fitted.values = Yhat, wi = wi,
        optim = fitCount, method = method, control = control, start = start,
        weights =
        if (identical(as.vector(weights), rep.int(1L, n))) NULL else weights,
        n = nobs, iterations = fitCount$niter, execTime = fitCount$xtimes, 
        loglik = fitCount$value, df.residual = nobs - length(coefs),
        vcov = varCovarcount, dist = dist,  link = linkList,
        converged = fitCount$convcode[1] == 0, data = data,
        formula = formula(f), call = cl, anc = anc, convPars = convPars,
        customPars = customPars, time = time, seriesPars = seriesPars,
        weiMethod = weiMethod
        )
    if (model) 
        rval$model <- mf
    if (y) 
        rval$y <- Y
    if (x) 
        rval$x <- modelMatrixList
    
    class(rval) <- "renewal"
    return(rval)
}
