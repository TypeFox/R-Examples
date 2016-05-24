predNames = c("algorithm", "score")
vbs <-
function(data=NULL) {
    if(is.null(data)) {
        stop("Need data to determine virtual best!")
    }
    lens = sapply(data$best, length)
    idxs = unlist(lapply(1:length(lens), function(i) { rep.int(i, lens[i]) }))
    ids = data$data[idxs,data$ids,drop=F]
    bests = unlist(data$best)
    scores = rep.int(1, length(bests))
    scores[is.na(bests)] = 0
    return(data.frame(ids, algorithm = bests, score = scores, iteration = 1))
}
class(vbs) = "llama.model"
attr(vbs, "type") = "virtual best"
attr(vbs, "hasPredictions") = FALSE
attr(vbs, "addCosts") = FALSE

singleBestByCount <-
function(data=NULL) {
    if(is.null(data)) {
        stop("Need data to determine single best!")
    }
    best = breakBestTies(data)
    ids = data$data[data$ids]
    data.frame(ids, algorithm = names(sort(table(best), decreasing=T)[1]), score = 1, iteration = 1)
}
class(singleBestByCount) = "llama.model"
attr(singleBestByCount, "type") = "single best"
attr(singleBestByCount, "hasPredictions") = FALSE
attr(singleBestByCount, "addCosts") = FALSE

singleBestByPar <-
function(data=NULL, factor=10) {
    if(is.null(data)) {
        stop("Need data to determine single best!")
    }
    if(is.null(data$success)) {
        stop("Need successes to compute PAR scores.")
    }

    best = names(sort(sapply(data$performance, function(x) {
        prs = data.frame(data$data[,data$ids,drop=F], algorithm = x, score = 1, iteration = 1)
        model = list(predictions=prs, ids=c("id"))
        class(model) = "llama.model"
        attr(model, "hasPredictions") = TRUE
        attr(model, "addCosts") = FALSE
        mean(parscores(data, model, factor=factor))
    })))[1]
    ids = data$data[data$ids]
    data.frame(ids, algorithm = best, score = 1, iteration = 1)
}
class(singleBestByPar) = "llama.model"
attr(singleBestByPar, "type") = "single best"
attr(singleBestByPar, "hasPredictions") = FALSE
attr(singleBestByPar, "addCosts") = FALSE

singleBestBySuccesses <-
function(data=NULL) {
    if(is.null(data)) {
        stop("Need data to determine single best!")
    }
    if(is.null(data$success)) {
        stop("Need successes to compute successes.")
    }

    best = names(sort(sapply(data$performance, function(x) {
        prs = data.frame(data$data[,data$ids,drop=F], algorithm = x, score = 1, iteration = 1)
        model = list(predictions=prs, ids=c("id"))
        class(model) = "llama.model"
        attr(model, "hasPredictions") = TRUE
        attr(model, "addCosts") = FALSE
        mean(successes(data, model))
    }), decreasing=T))[1]
    ids = data$data[data$ids]
    data.frame(ids, algorithm = best, score = 1, iteration = 1)
}
class(singleBestBySuccesses) = "llama.model"
attr(singleBestBySuccesses, "type") = "single best"
attr(singleBestBySuccesses, "hasPredictions") = FALSE
attr(singleBestBySuccesses, "addCosts") = FALSE

singleBest <-
function(data=NULL) {
    if(is.null(data)) {
        stop("Need data to determine single best!")
    }

    best = names(sort(sapply(data$performance, function(x) { mean(data$data[,x]) }), decreasing=!data$minimize))[1]
    ids = data$data[data$ids]
    data.frame(ids, algorithm = best, score = 1, iteration = 1)
}
class(singleBest) = "llama.model"
attr(singleBest, "type") = "single best"
attr(singleBest, "hasPredictions") = FALSE
attr(singleBest, "addCosts") = FALSE

breakBestTies <-
function(data=NULL, fold=NULL) {
    if(is.null(data)) {
        stop("Need data to break ties!")
    }

    if(is.null(fold)) {
        perfs = data$data[data$performance]
    } else {
        perfs = data$data[data$train[[fold]],][data$performance]
    }
    order = names(sort(sapply(perfs, mean), decreasing=!data$minimize))
    optfun = if(data$minimize) { which.min } else { which.max }
    best = factor(apply(perfs[order], 1, function(x) { order[optfun(unlist(x))] }))
    names(best) = NULL
    return(best)
}

predTable <-
function(predictions=NULL, bestOnly=TRUE) {
    if(is.null(predictions)) {
        stop("Need predictions to tabulate!")
    }
    if(bestOnly) {
        ids = setdiff(names(predictions), c(predNames))
        lvls = levels(predictions$algorithm)
        algorithms = factor(lvls[as.vector(by(predictions$algorithm, predictions[ids], head, 1))])
    } else {
        algorithms = predictions$algorithm
    }
    sort(table(algorithms), decreasing = TRUE)
}
