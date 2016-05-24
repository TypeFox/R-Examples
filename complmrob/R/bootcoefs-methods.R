#' Bootstrap the regression coefficients for a robust linear regression model
#'
#' This function provides an easy interface and useful output to bootstrapping the regression
#' coefficients of robust linear regression models
#'
#' The default method is to use fast and robust bootstrap as described in the paper by M. Salibian-Barrera, et al.
#' (see references). The other options are to bootstrap the residuals or to bootstrap cases (observations),
#' but the sampling distribution of the estimates from these methods can be numerically instable and take
#' longer to compute.
#'
#' @param object the model to bootstrap the coefficients from
#' @param R the number of bootstrap replicates.
#' @param method one of \code{"frb"} for fast and robust bootstrap, \code{"residuals"} to resample
#'      the residuals or \code{"cases"} to resample the cases.
#' @param ncpus the number of CPUs to utilize for bootstrapping.
#' @param cl a snow or parallel cluster to use for bootstrapping.
#' @param ... currently ignored.
#' @return A list of type \code{bootcoefs} for which \code{\link[=print.bootcoefs]{print}},
#'      \code{\link{summary}} and \code{\link[=plot.bootcoefs]{plot}} methods are available
#' @export
#' @references M. Salibian-Barrera, S. Aelst, and G. Willems. Fast and robust bootstrap. Statistical Methods and Applications, 17(1):41-71, 2008.
#' @examples
#' data <- data.frame(lifeExp = state.x77[, "Life Exp"], USArrests[ , -3])
#' mUSArr <- complmrob(lifeExp ~ ., data = data)
#' bc <- bootcoefs(mUSArr, R = 200) # the number of bootstrap replicates should
#'                                  # normally be higher!
#' summary(bc)
#' plot(bc) # for the model diagnostic plots
bootcoefs <- function(object, R = 999, method = c("frb", "residuals", "cases"), ncpus = NULL, cl = NULL, ...) {
    UseMethod("bootcoefs", object);
}

#' @describeIn bootcoefs For robust linear regression models with compositional data
#' @export
#' @importFrom boot boot
bootcoefs.complmrob <- function(object, R = 999, method = c("frb", "residuals", "cases"), ncpus = NULL, cl = NULL, ...) {
    #
    # Initialize auxiliary variables
    #
    method <- match.arg(method);

    clSetup <- setupCluster(ncpus, cl);

    bootParams <- list(
        R = quote(R),
        coefind = quote(object$coefind),
        parallel = quote(clSetup$parallel),
        ncpus = quote(length(clSetup$cl)),
        cl = quote(clSetup$cl)
    )

    bootres <- list();

    tryCatch({
        if(method == "frb") {
            bootres <- lapply(object$models, function(m, bootParams) {
                bootParams$data <- quote(model.frame(m));
                bootParams$control <- quote(bootStatFastControl(m));
                bootParams$statistic <- quote(bootStatFast);

                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                return(eval(bc));
            }, bootParams);

            if(object$intercept == TRUE) {
                m <- object$models[[1]];
                bootParams$coefind <- 1;
                bootParams$data <- quote(model.frame(m));
                bootParams$control <- quote(bootStatFastControl(m));
                bootParams$statistic <- quote(bootStatFast);

                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                bootres[["(Intercept)"]] <- eval(bc);
            }
        } else if(method == "residuals") {
            bootres <- lapply(object$models, function(m, bootParams) {
                respInd <- attr(m$terms, "response");

                bootParams$data <- quote(data.frame(model.frame(m)[ , -respInd, drop = FALSE],
                    fit = fitted(m), resid = residuals(m)));
                bootParams$weights = quote(m$rweights);
                bootParams$statistic <- quote(bootStatResiduals);
                bootParams$intercept <- quote(object$intercept);

                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                return(eval(bc));
            }, bootParams);

            if(object$intercept == TRUE) {
                m <- object$models[[1]];
                respInd <- attr(m$terms, "response");
                bootParams$coefind <- 1;
                bootParams$data <- quote(data.frame(model.frame(m)[ , -respInd, drop = FALSE],
                    fit = fitted(m), resid = residuals(m)));
                bootParams$weights = quote(m$rweights);
                bootParams$statistic <- quote(bootStatResiduals);
                bootParams$intercept <- quote(object$intercept);

                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                bootres[["(Intercept)"]] <- eval(bc);
            }
        } else {
            bootres <- lapply(object$models, function(m, bootParams) {
                bootParams$data <- quote(model.frame(m));
                bootParams$formula <- quote(formula(m$terms));
                bootParams$weights = quote(m$rweights);
                bootParams$statistic <- quote(bootStatCases);

                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                return(eval(bc));
            }, bootParams);

            if(object$intercept == TRUE) {
                m <- object$models[[1]];
                bootParams$coefind <- 1;
                bootParams$data <- quote(model.frame(m));
                bootParams$formula <- quote(formula(m$terms));
                bootParams$weights = quote(m$rweights);
                bootParams$statistic <- quote(bootStatCases);

                bc <- as.call(c(list(expression(boot::boot)[[1]]), bootParams))
                bootres[["(Intercept)"]] <- eval(bc);
            }
        }

        if(object$intercept == TRUE) {
            #reorder bootres so that the intercept is first
            bootres <- bootres[c(length(bootres), seq_len(length(bootres) - 1))];
        }
    }, error = function(e) {
        print(e);
    }, finally = {
        if(clSetup$needToShutdownCluster == TRUE) {
            parallel::stopCluster(clSetup$cl);
        }
    });

    ret <- list(
        bootres = bootres,
        model = object,
        R = R
    );

    class(ret) <- c("bootcoefs", "bccomplmrob");
    return(ret);
}

#' @describeIn bootcoefs For standard robust linear regression models
#' @export
#' @importFrom boot boot
bootcoefs.lmrob <- function(object, R = 999, method = c("frb", "residuals", "cases"), ncpus = NULL, cl = NULL, ...) {
    #
    # Initialize auxiliary variables
    #
    method <- match.arg(method);

    clSetup <- setupCluster(ncpus, cl);

    bootParams <- list(
        R = quote(R),
        coefind = quote(object$coefind),
        parallel = quote(clSetup$parallel),
        ncpus = quote(length(clSetup$cl)),
        cl = quote(clSetup$cl)
    )

    bootres = NULL;

    tryCatch({
        if(method == "frb") {
            bootres <- boot::boot(data = model.frame(object), statistic = bootStatFast,
                R = R, parallel = clSetup$parallel, ncpus = length(clSetup$cl), cl = clSetup$cl,
                coefind = seq_along(coef(object)), control = bootStatFastControl(object));
        } else if(method == "residuals") {
            respInd <- attr(object$terms, "response");
            tmpData <- data.frame(model.frame(object)[ , -respInd, drop = FALSE], fit = fitted(object), resid = residuals(object));

            bootres <- boot::boot(data = tmpData, statistic = bootStatResiduals, weights = object$rweights,
                R = R, parallel = clSetup$parallel, ncpus = length(clSetup$cl), cl = clSetup$cl,
                intercept = (attr(object$terms, "intercept") == 1), coefind = seq_along(coef(object)));
        } else {
            bootres <- boot::boot(data = model.frame(object), statistic = bootStatCases, weights = object$rweights,
                R = R, parallel = clSetup$parallel, ncpus = length(clSetup$cl), cl = clSetup$cl,
                coefind = seq_along(coef(object)), formula = formula(object$terms));
        }
    }, error = function(e) {
        print(e);
    }, finally = {
        if(clSetup$needToShutdownCluster == TRUE) {
            parallel::stopCluster(clSetup$cl);
        }
    });

    ret <- list(
        bootres = bootres,
        model = object,
        R = R
    );

    class(ret) <- c("bootcoefs", "bclmrob");
    return(ret);
}

#' @import parallel
setupCluster <- function(ncpus, cl) {
    ##
    ## Setup clusters (if any)
    ##
    needToShutdownCluster <- FALSE;
    parallel <- "no";
    if(!is.null(ncpus) && is.null(cl)) {
        if (.Platform$OS.type == "unix") {
            cl <- parallel::makeForkCluster(nnodes = ncpus);
        } else {
            cl <- parallel::makePSOCKcluster(names = ncpus);
        }
        needToShutdownCluster <- TRUE;
    }

    if(!is.null(cl)) {
        parallel <- "snow";
        tryCatch({
            parallel::clusterEvalQ(cl, {
                loadNamespace("robustbase");
            });
            parallel::clusterExport(cl, varlist = c("isomLR"));
        }, error = function(e) {
            if(needToShutdownCluster == TRUE) {
                parallel::stopCluster(cl);
                parallel <<- "no";
                needToShutdownCluster <<- FALSE;
            }
            cl <<- NULL;
        }, finally = function(...) {})
    }

    return(list(
        needToShutdownCluster = needToShutdownCluster,
        parallel = parallel,
        cl = cl
    ));
}
