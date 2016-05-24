#' Fit elastic net, LASSO or ridge regression model
#' 
#' Using the \pkg{glmnet} package implementation.
#' 
#' The \code{alpha} parameter of \code{\link[glmnet]{glmnet}} controls the type of
#' penalty. Use \code{0} (default) for lasso only, \code{1} for ridge only, or
#' an intermediate for a combination. This is typically the parameter to tune
#' on. The shrinkage, controlled by the \code{lambda} parameter, can be left
#' unspecified for internal tuning (works the same way as
#' \code{\link{fit_glmnet}}).
#' 
#' @param x Dataset.
#' @param y Response vector. Can be of many different types for solving
#'   different problems, see \code{\link[glmnet]{glmnet}}.
#' @param family Determines the the type of problem to solve. Auto detected if
#'   \code{y} is numeric or survival. See \code{\link{family}} for details.
#' @param nfolds See \code{\link[glmnet]{cv.glmnet}}.
#' @param foldid See \code{\link[glmnet]{cv.glmnet}}.
#' @param alpha Regularization parameter, see \code{\link[glmnet]{glmnet}}.
#' @param lambda Regularization parameter, see \code{\link[glmnet]{glmnet}}.
#' @param ... Sent to \code{\link{fit_glmnet}} or \code{\link[glmnet]{cv.glmnet}}.
#' @return Fitted elastic net model.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{predict_glmnet}},
#'   \code{\link{importance_glmnet}}, \code{\link{modeling_procedure}}
#' @export
fit_glmnet <- function(x, y, family, nfolds, foldid, alpha=1, lambda=NULL, ...){
    nice_require("glmnet", "is required to fit elastic net models")
    if(is.data.frame(x)){
        notify_once(id = "glmnet 'x' is not matrix",
                    "glmnet only takes data sets in matrix form. The conversion in the fitting function introduces an extra copy of the data set in the memory.",
                    fun = message)
        x <- as.matrix(x)
    }
    if(missing(family)){
        if(inherits(y, "Surv")) family <- "cox" else
        if(is.factor(y))
            if(length(levels(y)) == 2) family <- "binomial" else
            family <- "multinomial" else
        if(is.integer(y) & all(y >= 0)) family <- "poisson" else
        if(is.numeric(y)) family <- "gaussian" else
        stop("Could not auto detect glmnet family, see `?fit_glmnet`.")
    }
    if(family == "mgaussian")
        stop("The glmnet wrapper is not implemented for multivariate response so far.")
        #nice_require("survival")
    
    if(length(lambda) != 1){
        if(missing(nfolds)){
            nfolds <- if(missing(foldid)) 10 else max(foldid)
        }
        if(missing(foldid)){
            foldid <- apply(
                sapply(
                    resample("crossvalidation", y, nfold=nfolds, nrepeat=1),
                    as.integer
                ) == 0, 1, which)
            if(!is.integer(foldid))
                stop("Could not create `foldid`, check the contents of `y`.")
        }
    }

    details <- list(y.class = class(y),
                    y.levels = levels(y) # NULL if not applicable
    )
    if(length(lambda) != 1){
        # Tune lambda
        c(details,
          glmnet::cv.glmnet(x, y, family=family, nfolds=nfolds, foldid=foldid,
                            alpha=alpha, lambda=lambda, ...)
        )
        #))
    } else {
        # Train single glmnet
        c(details,
          glmnet.fit = glmnet::glmnet(x, y, family=family,
                                      alpha=alpha, lambda=lambda, ...),
          list(lambda = lambda, lambda.min = lambda)
        )
    }
}


#' Predict using generalized linear model with elastic net regularization
#'
#' Due to the way \code{\link[glmnet]{glmnet}} is implemented, the regularization alpha
#' can not be modified after the model is fitted.
#'
#' @param object Fitted model.
#' @param x New data to be predicted.
#' @param s Regularization parameter lambda.
#' @param ... Sent to \code{\link[glmnet]{predict.glmnet}}.
#' @return A list with a subset of the following elements:
#' \describe{
#'     \item{\code{prediction}}{The response of the modeling problem, i.e. a
#'         factor for classification, problems, a numeric for regressions, and a
#'         relative risk for survival analyses.}
#'     \item{\code{probability}}{Data frame of predicted class probabilities.}
#'     \item{\code{link}}{Link function values.}
#' }
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{fit_glmnet}},
#'   \code{\link{importance_glmnet}}, \code{\link{modeling_procedure}}
#' @export
predict_glmnet <- function(object, x, s, ...){
    nice_require("glmnet", "is required to make precdictions with an elastic net model")
    if(missing(s)){
        if("lambda.min" %in% names(object)){
            s <- object$lambda.min
        } else {
            stop("The glmnet model has not been tuned.")
        }
    }
    predictor <- function(type=c("link", "class", "response"), ...){
        p <- predict(object$glmnet.fit, x, s=s, type=match.arg(type), ...)
        if(ncol(p) == 1) as.vector(p)
        else as.data.frame(p)
    }
    if(inherits(object$glmnet.fit, c("lognet", "multnet"))){
        # Classification
        p <- list(prediction = predictor("class", ...),
                  probability = predictor("response", ...))
        if(inherits(object$glmnet.fit, "lognet")){
            p$probability <- data.frame(V1=1-p$probability,
                                        V2=p$probability)
        }
        if(!is.null(object$y.levels)){
            p$prediction %<>% factor(object$y.levels)
            colnames(p$probability) <- object$y.levels
        }
        p
    } else if(inherits(object$glmnet.fit, "elnet")){ 
        # Regression
        list(prediction = predictor())
    } else if(inherits(object$glmnet.fit, "mrelnet")){
        stop("The glmnet wrapper is not implemented for multivariate response so far.")
    } else if(inherits(object$glmnet.fit, c("coxnet", "fishnet"))){
        # Cox regression (for survival analysis)
        # Non-negative counts regression
        list(prediction = predictor("response"),
             link = predictor("link"))
    }
}

#' Feature importance extractor for elastic net models
#'
#' @param object Fitted elastic net model, as produced by
#'   \code{\link{fit_glmnet}}.
#' @param s Regularization parameter lambda.
#' @param ... Sent to \code{\link[glmnet]{predict.glmnet}}.
#' @return A feature imortance data frame.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{fit_glmnet}},
#'   \code{\link{predict_glmnet}}, \code{\link{modeling_procedure}}
#' @export
importance_glmnet <- function(object, s, ...){
    if(missing(s)){
        if("lambda.min" %in% names(object)){
            s <- object$lambda.min
        } else {
            stop("The glmnet model has not been tuned.")
        }
    }
    imp <- predict(object$glmnet.fit, s=s, type="coef", ...)
    if(is.list(imp))
        imp <- do.call(cbind, imp)
    imp <- as.matrix(imp[-1,])
    colnames(imp) <- if(!is.null(object$y.level)){
        object$y.level
    } else if (ncol(imp) == 1){
        "coefficient"
    } else {
        paste("coefficient", 1:ncol(imp), sep="")
    }
    data.frame(feature = rownames(imp), imp, row.names=NULL)
}

#' @rdname fit_glmnet
#' @export
fit_ridge_regression <- function(...){
    fit_glmnet(..., alpha=0)
}
#' @rdname predict_glmnet
#' @export
predict_ridge_regression <- predict_glmnet
#' @rdname importance_glmnet
#' @export
importance_ridge_regression <- importance_glmnet

#' @rdname fit_glmnet
#' @export
fit_lasso <- function(...){
    fit_glmnet(..., alpha=1)
}
#' @rdname predict_glmnet
#' @export
predict_lasso <- predict_glmnet
#' @rdname importance_glmnet
#' @export
importance_lasso <- importance_glmnet

