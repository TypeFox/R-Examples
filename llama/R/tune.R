# code adapted from aslib package
tuneModel <-
function(ldf, llama.fun, learner, design, metric = parscores, nfolds = 10L, quiet = FALSE) {
    assertClass(ldf, "llama.data")
    assertClass(llama.fun, "llama.modelFunction")
    assertClass(learner, "Learner")
    assertClass(design, "data.frame")
    assertClass(metric, "llama.metric")
    assertInteger(nfolds)

    if(length(attr(ldf, "hasSplits")) == 0) {
        n.outer.folds = nfolds
        ldf = cvFolds(ldf, nfolds = nfolds)
    } else {
        n.outer.folds = length(ldf$test)
    }
    inner.retval = lapply(1:n.outer.folds, function(i) {
        if(!quiet) message(paste("Fold ", i, "/", n.outer.folds, ":", sep = ""))
        ldf2 = ldf
        ldf2$data = ldf$data[ldf$train[[i]],]
        ldf3 = ldf2
        ldf3 = cvFolds(ldf3, nfolds = nfolds)
        best.parvals = tuneLlamaModel(ldf3, llama.fun, learner, design, metric, quiet)

        outer.split.ldf = ldf
        outer.split.ldf$train = list(ldf$train[[i]])
        outer.split.ldf$test = list(ldf$test[[i]])

        learner2 = setHyperPars(learner, par.vals = best.parvals)
        model = llama.fun(learner2, data = outer.split.ldf)

        retval = model$predictions
        retval$iteration = i
        return(list(predictions = retval, parvals = best.parvals))
    })

    best.parvals = tuneLlamaModel(ldf, llama.fun, learner, design, metric, quiet)
    learner2 = setHyperPars(learner, par.vals = best.parvals)
    full.split.ldf = ldf
    full.split.ldf$train = list(ldf$train[[1]])
    full.split.ldf$test = list(ldf$test[[1]])
    model = llama.fun(learner2, data = full.split.ldf)

    predictions = rbind.fill(lapply(inner.retval, function(x) x$predictions))
    parvals = lapply(inner.retval, function(x) x$parvals)

    retval = list(predictions = predictions, models = model$models, predictor = model$predictor,
        parvals = best.parvals, inner.parvals = parvals)
    class(retval) = "llama.model"
    attr(retval, "type") = attr(model, "type")
    attr(retval, "hasPredictions") = TRUE
    attr(retval, "addCosts") = TRUE

    return(retval)
}

tuneLlamaModel <-
function(ldf, llama.fun, learner, design, metric, quiet) {
    # FIXME: we currently do not handle failed tuning evals
    ys = parallelMap(function(x) {
        pars = as.list(design[x,,drop = FALSE])
        learner = setHyperPars(learner, par.vals = pars)
        model = llama.fun(learner, ldf)
        score = mean(metric(ldf, model))

        if(!quiet) message(paste("      [", paste(names(pars), pars, sep = " = ", collapse = ", "), "], score = ", score, sep = ""))
        return(score)
    }, 1:nrow(design), simplify = TRUE, level = "llama.tune")

    if(attr(metric, "minimize")) {
        best.i = getMinIndex(ys)
    } else {
        best.i = getMaxIndex(ys)
    }
    best.parvals = as.list(design[best.i,,drop = FALSE])
    if(!quiet) message(paste("Best: [", paste(names(best.parvals), best.parvals, sep = " = ", collapse = ", "), "], score = ", ys[best.i], sep = ""))

    return(best.parvals)
}
