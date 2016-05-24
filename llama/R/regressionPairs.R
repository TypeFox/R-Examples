regressionPairs <-
function(regressor=NULL, data=NULL, pre=function(x, y=NULL) { list(features=x) }, combine=NULL, save.models=NA, use.weights = TRUE) {
    if(!testClass(regressor, "Learner")) {
        stop("Need regressor!")
    }
    assertClass(data, "llama.data")
    hs = attr(data, "hasSplits")
    if(is.null(hs) || hs != TRUE) {
        stop("Need data with train/test split!")
    }

    signPair = function(k, p, sign, row, combns) {
        if(p == combns[1,k]) {
            sign(row[k])
        } else if(p == combns[2,k]) {
            sign(-row[k])
        } else {
            0
        }
    }

    worstScore = if(data$minimize) { Inf } else { -Inf }

    totalBests = data.frame(target=factor(breakBestTies(data), levels=data$performance))
    combns = combn(data$performance, 2)
    predictions = rbind.fill(parallelMap(function(i) {
        trf = pre(data$data[data$train[[i]],][data$features])
        tsf = pre(data$data[data$test[[i]],][data$features], trf$meta)
        ids = data$data[data$test[[i]],][data$ids]
        trp = data$data[data$train[[i]],][data$performance]

        trainpredictions = matrix(nrow=nrow(trf$features), ncol=ncol(combns))
        pairpredictions = matrix(nrow=nrow(tsf$features), ncol=ncol(combns))
        for (j in 1:ncol(combns)) {
            values = data.frame(target=data$data[data$train[[i]],combns[1,j]] - data$data[data$train[[i]],combns[2,j]])
            task = makeRegrTask(id="regression", target="target", data=cbind(values, trf$features), fixup.data="quiet")
            model = train(regressor, task = task)
            if(!is.na(save.models)) {
                saveRDS(list(model=model, train.data=task, test.data=tsf$features), file = paste(save.models, regressor$id, combns[1,j], combns[2,j], i, "rds", sep="."))
            }
            if(!is.null(combine)) { # only do this if we need it
                trainpredictions[,j] = getPredictionResponse(predict(model, newdata=trf$features))
            }
            pairpredictions[,j] = getPredictionResponse(predict(model, newdata=tsf$features))
        }

        if(!is.null(combine)) {
            trainBests = data.frame(target=factor(breakBestTies(data, i), levels=data$performance))
            if(hasLearnerProperties(combine, "weights") && use.weights) {
                trw = abs(apply(trp, 1, max) - apply(trp, 1, min))
                task = makeClassifTask(id="regression", target="target", weights=trw, data=cbind(trainBests, trf$features, data.frame(trainpredictions)), fixup.data="quiet", check.data=FALSE)
            } else {
                task = makeClassifTask(id="regression", target="target", data=cbind(trainBests, trf$features, data.frame(trainpredictions)), fixup.data="quiet", check.data=FALSE)
            }
            combinedmodel = train(combine, task = task)
            if(!is.na(save.models)) {
                saveRDS(list(model=combinedmodel, train.data=task, test.data=cbind(tsf$features, data.frame(pairpredictions))), file = paste(save.models, combine$id, "combined", i, "rds", sep="."))
            }
            preds = getPredictionResponse(predict(combinedmodel, newdata=cbind(tsf$features, data.frame(pairpredictions))))
            combinedpredictions = rbind.fill(lapply(1:length(preds), function(j) {
                if(all(is.na(preds[j,drop=F]))) {
                    data.frame(ids[j,,drop=F], algorithm=NA, score=worstScore, iteration=i, row.names = NULL)
                } else {
                    tab = table(preds[j,drop=F])
                    data.frame(ids[j,,drop=F], algorithm=names(tab), score=as.vector(tab), iteration=i, row.names = NULL)
                }
            }))
        } else {
            if(data$minimize) {
                sign = function(d) { return(-d) }
            } else {
                sign = function(d) { return(d) }
            }
            combinedpredictions = rbind.fill(lapply(1:nrow(pairpredictions), function(j) {
                row = pairpredictions[j,]
                performanceSums = sapply(data$performance, function(p) {
                    sum(sapply(1:ncol(combns), signPair, p, sign, row, combns))
                })
                x = sort(performanceSums, decreasing = TRUE)
                if(all(is.na(x))) {
                    data.frame(ids[j,,drop=F], algorithm=NA, score=worstScore, iteration=i, row.names = NULL)
                } else {
                    data.frame(ids[j,,drop=F], algorithm=names(x), score=unlist(x), iteration=i, row.names = NULL)
                }
            }))
        }
        return(combinedpredictions)
    }, 1:length(data$train), level = "llama.fold"))

    fs = pre(data$data[data$features])
    fp = data$data[data$performance]
    fw = abs(apply(fp, 1, max) - apply(fp, 1, min))
    models = lapply(1:ncol(combns), function(i) {
        values = data.frame(target=data$data[[combns[1,i]]] - data$data[[combns[2,i]]])
        task = makeRegrTask(id="regression", target="target", data=cbind(values, fs$features), fixup.data="quiet")
        return(train(regressor, task = task))
    })
    if(!is.null(combine)) {
        trainpredictions = matrix(nrow=nrow(fs$features), ncol=ncol(combns))
        for(i in 1:ncol(combns)) {
            trainpredictions[,i] = getPredictionResponse(predict(models[[i]], newdata=fs$features))
        }
        if(hasLearnerProperties(combine, "weights") && use.weights) {
            task = makeClassifTask(id="regression", target="target", weights=fw, data=cbind(totalBests, fs$features, data.frame(trainpredictions)), fixup.data="quiet", check.data=FALSE)
        } else {
            task = makeClassifTask(id="regression", target="target", data=cbind(totalBests, fs$features, data.frame(trainpredictions)), fixup.data="quiet", check.data=FALSE)
        }
        combinedmodel = train(combine, task = task)
    }

    predictor = function(x) {
        tsf = pre(x[data$features], fs$meta)
        if(length(intersect(colnames(x), data$ids)) > 0) {
            ids = x[data$ids]
        } else {
            ids = data.frame(id = 1:nrow(x)) # don't have IDs, generate them
        }
        pairpredictions = matrix(nrow=nrow(tsf$features), ncol=ncol(combns))
        for(i in 1:ncol(combns)) {
            pairpredictions[,i] = getPredictionResponse(predict(models[[i]], newdata=tsf$features))
        }
        if(!is.null(combine)) {
            preds = getPredictionResponse(predict(combinedmodel, newdata=cbind(tsf$features, data.frame(pairpredictions))))
            combinedpredictions = rbind.fill(lapply(1:length(preds), function(j) {
                if(all(is.na(preds[j,drop=F]))) {
                    data.frame(ids[j,,drop=F], algorithm=NA, score=worstScore, iteration=i, row.names = NULL)
                } else {
                    tab = table(preds[j,drop=F])
                    data.frame(ids[j,,drop=F], algorithm=names(tab), score=as.vector(tab), iteration=1, row.names = NULL)
                }
            }))
        } else {
            if(data$minimize) {
                sign = function(d) { return(-d) }
            } else {
                sign = function(d) { return(d) }
            }
            combinedpredictions = rbind.fill(lapply(1:nrow(pairpredictions), function(j) {
                row = pairpredictions[j,]
                performanceSums = sapply(data$performance, function(p) {
                    sum(sapply(1:ncol(combns), signPair, p, sign, row, combns))
                })
                x = sort(performanceSums, decreasing = TRUE)
                if(all(is.na(x))) {
                    data.frame(ids[j,,drop=F], algorithm=NA, score=worstScore, iteration=1, row.names = NULL)
                } else {
                    data.frame(ids[j,,drop=F], algorithm=names(x), score=unlist(x), iteration=1, row.names = NULL)
                }
            }))
        }
        return(combinedpredictions)
    }
    class(predictor) = "llama.model"
    attr(predictor, "type") = "regressionPairs"
    attr(predictor, "hasPredictions") = FALSE
    attr(predictor, "addCosts") = TRUE

    retval = list(predictions=predictions, models=models, predictor=predictor)
    class(retval) = "llama.model"
    attr(retval, "type") = "regressionPairs"
    attr(retval, "hasPredictions") = TRUE
    attr(retval, "addCosts") = TRUE

    return(retval)
}
class(regressionPairs) = "llama.modelFunction"
