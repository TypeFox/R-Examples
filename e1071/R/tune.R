tune.control <- function(random = FALSE,
                         nrepeat = 1,
                         repeat.aggregate = mean,
                         sampling = c("cross", "fix", "bootstrap"),
                         sampling.aggregate = mean,
                         sampling.dispersion = sd,
                         cross = 10,
                         fix = 2 / 3,
                         nboot = 10,
                         boot.size = 9 / 10,
                         best.model = TRUE,
                         performances = TRUE,
                         error.fun = NULL) {
    structure(list(random = random,
                   nrepeat = nrepeat,
                   repeat.aggregate = repeat.aggregate,
                   sampling = match.arg(sampling),
                   sampling.aggregate = sampling.aggregate,
                   sampling.dispersion = sampling.dispersion,
                   cross = cross,
                   fix = fix,
                   nboot = nboot,
                   boot.size = boot.size,
                   best.model = best.model,
                   performances = performances,
		   error.fun = error.fun
                   ),
              class = "tune.control"
              )
}

tune <- function(method, train.x, train.y = NULL, data = list(),
                 validation.x = NULL, validation.y = NULL,
                 ranges = NULL, predict.func = predict,
                 tunecontrol = tune.control(),
                 ...
                 ) {
    call <- match.call()

    ## internal helper functions
    resp <- function(formula, data) {

        model.response(model.frame(formula, data))
    }

    classAgreement <- function (tab) {
        n <- sum(tab)
        if (!is.null(dimnames(tab))) {
            lev <- intersect(colnames(tab), rownames(tab))
            p0 <- sum(diag(tab[lev, lev])) / n
        } else {
            m <- min(dim(tab))
            p0 <- sum(diag(tab[1:m, 1:m])) / n
        }
        p0
    }

    ## parameter handling
    if (tunecontrol$sampling == "cross")
        validation.x <- validation.y <- NULL
    useFormula <- is.null(train.y)
    if (useFormula && (is.null(data) || length(data) == 0))
        data <- model.frame(train.x)
    if (is.vector(train.x)) train.x <- t(t(train.x))
    if (is.data.frame(train.y))
        train.y <- as.matrix(train.y)

    ## prepare training indices
    if (!is.null(validation.x)) tunecontrol$fix <- 1
    n <- nrow(if (useFormula) data else train.x)
    perm.ind <- sample(n)
    if (tunecontrol$sampling == "cross") {
        if (tunecontrol$cross > n)
            stop(sQuote("cross"), " must not exceed sampling size!")
        if (tunecontrol$cross == 1)
            stop(sQuote("cross"), " must be greater than 1!")
    }
    train.ind <- if (tunecontrol$sampling == "cross")
        tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
    else if (tunecontrol$sampling == "fix")
        list(perm.ind[1:trunc(n * tunecontrol$fix)])
    else ## bootstrap
        lapply(1:tunecontrol$nboot,
               function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))

    ## find best model
    parameters <- if (is.null(ranges))
        data.frame(dummyparameter = 0)
    else
        expand.grid(ranges)
    p <- nrow(parameters)
    if (!is.logical(tunecontrol$random)) {
        if (tunecontrol$random < 1)
            stop("random must be a strictly positive integer")
        if (tunecontrol$random > p) tunecontrol$random <- p
        parameters <- parameters[sample(1:p, tunecontrol$random),]
    }
    model.variances <- model.errors <- c()

    ## - loop over all models
    for (para.set in 1:p) {
        sampling.errors <- c()

        ## - loop over all training samples
        for (sample in 1:length(train.ind)) {
            repeat.errors <- c()

            ## - repeat training `nrepeat' times
            for (reps in 1:tunecontrol$nrepeat) {

                ## train one model
                pars <- if (is.null(ranges))
                    NULL
                else
                    lapply(parameters[para.set,,drop = FALSE], unlist)

                model <- if (useFormula)
                    do.call(method, c(list(train.x,
                                           data = data,
                                           subset = train.ind[[sample]]),
                                      pars, list(...)
                                      )
                            )
                else
                    do.call(method, c(list(train.x[train.ind[[sample]],],
                                           y = train.y[train.ind[[sample]]]),
                                      pars, list(...)
                                      )
                            )

                ## predict validation set
                pred <- predict.func(model,
                                     if (!is.null(validation.x))
                                     validation.x
                                     else if (useFormula)
                                     data[-train.ind[[sample]],,drop = FALSE]
                                     else if (inherits(train.x, "matrix.csr"))
                                     train.x[-train.ind[[sample]],]
                                     else
                                     train.x[-train.ind[[sample]],,drop = FALSE]
                                     )

                ## compute performance measure
                true.y <- if (!is.null(validation.y))
                    validation.y
                else if (useFormula) {
                    if (!is.null(validation.x))
                        resp(train.x, validation.x)
                    else
                        resp(train.x, data[-train.ind[[sample]],])
                } else
                    train.y[-train.ind[[sample]]]

                if (is.null(true.y)) true.y <- rep(TRUE, length(pred))

                repeat.errors[reps] <- if (!is.null(tunecontrol$error.fun))
                    tunecontrol$error.fun(true.y, pred)
                else if ((is.logical(true.y) || is.factor(true.y)) && (is.logical(pred) || is.factor(pred) || is.character(pred))) ## classification error
                    1 - classAgreement(table(pred, true.y))
                else if (is.numeric(true.y) && is.numeric(pred)) ## mean squared error
                    crossprod(pred - true.y) / length(pred)
                else
                    stop("Dependent variable has wrong type!")
            }
            sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)
        }
        model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors)
        model.variances[para.set] <- tunecontrol$sampling.dispersion(sampling.errors)
    }

    ## return results
    best <- which.min(model.errors)
    pars <- if (is.null(ranges))
        NULL
    else
        lapply(parameters[best,,drop = FALSE], unlist)
    structure(list(best.parameters  = parameters[best,,drop = FALSE],
                   best.performance = model.errors[best],
                   method           = if (!is.character(method))
                   deparse(substitute(method)) else method,
                   nparcomb         = nrow(parameters),
                   train.ind        = train.ind,
                   sampling         = switch(tunecontrol$sampling,
                   fix = "fixed training/validation set",
                   bootstrap = "bootstrapping",
                   cross = if (tunecontrol$cross == n) "leave-one-out" else
                   paste(tunecontrol$cross,"-fold cross validation", sep="")
                   ),
                   performances     = if (tunecontrol$performances) cbind(parameters, error = model.errors, dispersion = model.variances),
                   best.model       = if (tunecontrol$best.model) {
                       modeltmp <- if (useFormula)
                           do.call(method, c(list(train.x, data = data),
                                             pars, list(...)))
                       else
                           do.call(method, c(list(x = train.x,
                                                  y = train.y),
                                             pars, list(...)))
                       call[[1]] <- as.symbol("best.tune")
                       modeltmp$call <- call
                       modeltmp
                   }
                   ),
              class = "tune"
              )
}

best.tune <- function(...) {
    call <- match.call()
    modeltmp <- tune(...)$best.model
    modeltmp$call <- call
    modeltmp
}

print.tune <- function(x, ...) {
    if (x$nparcomb > 1) {
        cat("\nParameter tuning of ", sQuote(x$method), ":\n\n", sep="")
        cat("- sampling method:", x$sampling,"\n\n")
        cat("- best parameters:\n")
        tmp <- x$best.parameters
        rownames(tmp) <- ""
        print(tmp)
        cat("\n- best performance:", x$best.performance, "\n")
        cat("\n")
    } else {
        cat("\nError estimation of ", sQuote(x$method), " using ", x$sampling, ": ",
            x$best.performance, "\n\n", sep="")
    }
}

summary.tune <- function(object, ...)
    structure(object, class = "summary.tune")

print.summary.tune <- function(x, ...) {
    print.tune(x)
    if (!is.null(x$performances) && (x$nparcomb > 1)) {
        cat("- Detailed performance results:\n")
        print(x$performances)
        cat("\n")
    }
}

hsv_palette <- function(h = 2/3, from = 0.7, to = 0.2, v = 1)
    function(n) hsv(h = h, s = seq(from, to, length.out = n), v = v)

plot.tune <- function(x,
                      type=c("contour","perspective"),
                      theta=60,
                      col="lightblue",
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      swapxy = FALSE,
                      transform.x = NULL,
                      transform.y = NULL,
                      transform.z = NULL,
                      color.palette = hsv_palette(),
                      nlevels = 20,
                      ...)
{
    if (is.null(x$performances))
        stop("Object does not contain detailed performance measures!")
    k <- ncol(x$performances)
    if (k > 4) stop("Cannot visualize more than 2 parameters")
    type = match.arg(type)

    if (is.null(main))
        main <- paste("Performance of `", x$method, "'", sep="")

    if (k == 3)
        plot(x$performances[,1:2], type = "b", main = main)
    else  {
        if (!is.null(transform.x))
            x$performances[,1] <- transform.x(x$performances[,1])
        if (!is.null(transform.y))
            x$performances[,2] <- transform.y(x$performances[,2])
        if (!is.null(transform.z))
            x$performances[,3] <- transform.z(x$performances[,3])
        if (swapxy)
            x$performances[,1:2] <- x$performances[,2:1]
        x <- xtabs(error~., data = x$performances[,-k])
        if (is.null(xlab)) xlab <- names(dimnames(x))[1 + swapxy]
        if (is.null(ylab)) ylab <- names(dimnames(x))[2 - swapxy]
        if (type == "perspective")
            persp(x=as.double(rownames(x)),
                  y=as.double(colnames(x)),
                  z=x,
                  xlab=xlab,
                  ylab=ylab,
                  zlab="accuracy",
                  theta=theta,
                  col=col,
                  ticktype="detailed",
                  main = main,
                  ...
                  )
        else
            filled.contour(x=as.double(rownames(x)),
                           y=as.double(colnames(x)),
                           xlab=xlab,
                           ylab=ylab,
                           nlevels=nlevels,
                           color.palette = color.palette,
                           main = main,
                           x, ...)
    }
}

#############################################
## convenience functions for some methods
#############################################

tune.svm <- function(x, y = NULL, data = NULL, degree = NULL, gamma = NULL,
                     coef0 = NULL, cost = NULL, nu = NULL, class.weights = NULL,
                     epsilon = NULL, ...) {
    call <- match.call()
    call[[1]] <- as.symbol("best.svm")
    ranges <- list(degree = degree, gamma = gamma,
                   coef0 = coef0, cost = cost, nu = nu,
                   class.weights = class.weights, epsilon = epsilon)
    ranges[sapply(ranges, is.null)] <- NULL
    if (length(ranges) < 1)
        ranges = NULL
    modeltmp <- if (inherits(x, "formula"))
        tune("svm", train.x = x, data = data, ranges = ranges, ...)
    else
        tune("svm", train.x = x, train.y = y, ranges = ranges, ...)
    if (!is.null(modeltmp$best.model))
        modeltmp$best.model$call <- call
    modeltmp
}

best.svm <- function(x, tunecontrol = tune.control(), ...) {
    call <- match.call()
    tunecontrol$best.model = TRUE
    modeltmp <- tune.svm(x, ..., tunecontrol = tunecontrol)$best.model
    modeltmp$call <- call
    modeltmp
}

tune.nnet <- function(x, y = NULL, data = NULL,
                      size = NULL, decay = NULL, trace = FALSE,
                      tunecontrol = tune.control(nrepeat = 5),
                      ...) {
    call <- match.call()
    call[[1]] <- as.symbol("best.nnet")
    loadNamespace("nnet")
    predict.func <- predict
    useFormula <- inherits(x, "formula")
    if (is.factor(y) ||
        (useFormula && is.factor(model.response(model.frame(formula = x, data = data))))
        )
        predict.func = function(...) predict(..., type = "class")
    ranges <- list(size = size, decay = decay)
    ranges[sapply(ranges, is.null)] <- NULL
    if (length(ranges) < 1)
        ranges = NULL
    modeltmp <- if (useFormula)
        tune("nnet", train.x = x, data = data, ranges = ranges, predict.func = predict.func,
             tunecontrol = tunecontrol, trace = trace, ...)
    else
        tune("nnet", train.x = x, train.y = y, ranges = ranges, predict.func = predict.func,
             tunecontrol = tunecontrol, trace = trace, ...)
    if (!is.null(modeltmp$best.model))
        modeltmp$best.model$call <- call
    modeltmp
}

best.nnet <- function(x, tunecontrol = tune.control(nrepeat = 5), ...) {
    call <- match.call()
    tunecontrol$best.model = TRUE
    modeltmp <- tune.nnet(x, ..., tunecontrol = tunecontrol)$best.model
    modeltmp$call <- call
    modeltmp
}

tune.randomForest <- function(x, y = NULL, data = NULL, nodesize = NULL, mtry = NULL, ntree = NULL, ...) {
    call <- match.call()
    call[[1]] <- as.symbol("best.randomForest")
    loadNamespace("randomForest")
    ranges <- list(nodesize = nodesize, mtry = mtry, ntree = ntree)
    ranges[sapply(ranges, is.null)] <- NULL
    if (length(ranges) < 1)
        ranges = NULL
    modeltmp <- if (inherits(x, "formula"))
        tune("randomForest", train.x = x, data = data, ranges = ranges, ...)
    else
        tune("randomForest", train.x = x, train.y = y, ranges = ranges, ...)
    if (!is.null(modeltmp$best.model))
        modeltmp$best.model$call <- call
    modeltmp
}

best.randomForest <- function(x, tunecontrol = tune.control(), ...) {
    call <- match.call()
    tunecontrol$best.model = TRUE
    modeltmp <- tune.randomForest(x, ..., tunecontrol = tunecontrol)$best.model
    modeltmp$call <- call
    modeltmp
}

knn.wrapper <- function(x, y, k = 1, l = 0, ...)
    list(train = x, cl = y, k = k, l = l, ...)

tune.knn <- function(x, y, k = NULL, l = NULL, ...) {
    loadNamespace("class")
    ranges <- list(k = k, l = l)
    ranges[sapply(ranges, is.null)] <- NULL
    if (length(ranges) < 1)
        ranges = NULL
    tune("knn.wrapper",
         train.x = x, train.y = y, ranges = ranges,
         predict.func = function(x, ...) knn(train = x$train, cl = x$cl, k = x$k, l = x$l, ...),
         ...)
}

rpart.wrapper <- function(formula, minsplit=20, minbucket=round(minsplit/3), cp=0.01,
                          maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10,
                          surrogatestyle=0, maxdepth=30, ...)
    rpart::rpart(formula,
                 control = rpart::rpart.control(minsplit=minsplit, minbucket=minbucket,
                 cp=cp, maxcompete=maxcompete, maxsurrogate=maxsurrogate,
                 usesurrogate=usesurrogate, xval=xval,
                 surrogatestyle=surrogatestyle, maxdepth=maxdepth),
                 ...
                 )

tune.rpart <- function(formula, data, na.action = na.omit,
                       minsplit=NULL, minbucket=NULL, cp=NULL,
                       maxcompete=NULL, maxsurrogate=NULL, usesurrogate=NULL, xval=NULL,
                       surrogatestyle=NULL, maxdepth=NULL,
                       predict.func = NULL,
                       ...) {
    call <- match.call()
    call[[1]] <- as.symbol("best.rpart")
    loadNamespace("rpart")
    ranges <- list(minsplit=minsplit, minbucket=minbucket, cp=cp,
                   maxcompete=maxcompete, maxsurrogate=maxsurrogate,
                   usesurrogate=usesurrogate, xval=xval,
                   surrogatestyle=surrogatestyle, maxdepth=maxdepth)
    ranges[sapply(ranges, is.null)] <- NULL
    if (length(ranges) < 1)
        ranges <- NULL

    predict.func <- if (is.factor(model.response(model.frame(formula, data))))
        function(...) predict(..., type = "class")
    else
        predict
    modeltmp <- tune("rpart.wrapper", train.x = formula, data = data, ranges = ranges,
                     predict.func = predict.func, na.action = na.action, ...)
    if (!is.null(modeltmp$best.model))
        modeltmp$best.model$call <- call
    modeltmp
}

best.rpart <- function(formula, tunecontrol = tune.control(), ...) {
    call <- match.call()
    tunecontrol$best.model = TRUE
    modeltmp <- tune.rpart(formula, ..., tunecontrol = tunecontrol)$best.model
    modeltmp$call <- call
    modeltmp
}
