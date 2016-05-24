classifyPairs <-
function(classifier=NULL, data=NULL, pre=function(x, y=NULL) { list(features=x) }, combine=NULL, save.models=NA, use.weights = TRUE) {
    if(!testClass(classifier, "Learner")) {
        stop("Need classifier!")
    }
    assertClass(data, "llama.data")
    hs = attr(data, "hasSplits")
    if(is.null(hs) || hs != TRUE) {
        stop("Need data with train/test split!")
    }

    totalBests = data.frame(target=factor(breakBestTies(data), levels=data$performance))
    combns = combn(data$performance, 2)
    predictions = rbind.fill(parallelMap(function(i) {
        trf = pre(data$data[data$train[[i]],][data$features])
        tsf = pre(data$data[data$test[[i]],][data$features], trf$meta)
        ids = data$data[data$test[[i]],][data$ids]
        trp = data$data[data$train[[i]],][data$performance]

        trainpredictions = list()
        pairpredictions = list()
        for (j in 1:ncol(combns)) {
            if(data$minimize) {
                cmp = function(x, y) {
                    sapply(data$data[data$train[[i]],][x] < data$data[data$train[[i]],][y], function(z) { if(z) { x } else { y } })
                }
            } else {
                cmp = function(x, y) {
                    sapply(data$data[data$train[[i]],][x] > data$data[data$train[[i]],][y], function(z) { if(z) { x } else { y } })
                }
            }
            labels = data.frame(target=factor(cmp(combns[1,j], combns[2,j])))
            if(hasLearnerProperties(classifier, "weights") && use.weights) {
                trw = abs(data$data[data$train[[i]],combns[1,j]] - data$data[data$train[[i]],combns[2,j]])
                task = makeClassifTask(id="classifyPairs", target="target", weights=trw, data=data.frame(labels, trf$features), fixup.data="quiet", check.data=FALSE)
            } else {
                task = makeClassifTask(id="classifyPairs", target="target", data=data.frame(labels, trf$features), fixup.data="quiet", check.data=FALSE)
            }
            if(length(unique(labels$target)) == 1) {
                # one-class problem
                model = train(constantClassifier, task = task)
            } else {
                model = train(classifier, task = task)
            }
            if(!is.na(save.models)) {
                saveRDS(list(model=model, train.data=task, test.data=tsf$features), file = paste(save.models, classifier$id, combns[1,j], combns[2,j], i, "rds", sep="."))
            }
            if(!is.null(combine)) { # only do this if we need it
                preds = predict(model, newdata=trf$features)
                trainpredictions[[j]] = if(preds$predict.type == "prob") {
                    getPredictionProbabilities(preds, preds$task.desc$class.levels)
                } else {
                    tmp = getPredictionResponse(preds)
                    rbind.fill(lapply(tmp, function(x) data.frame(t(setNames(as.numeric(x == levels(tmp)), levels(tmp))))))
                }
            }
            preds = predict(model, newdata=tsf$features)
            pairpredictions[[j]] = if(preds$predict.type == "prob") {
                getPredictionProbabilities(preds, preds$task.desc$class.levels)
            } else {
                tmp = getPredictionResponse(preds)
                rbind.fill(lapply(tmp, function(x) data.frame(t(setNames(as.numeric(x == levels(tmp)), levels(tmp))))))
            }
        }

        if(!is.null(combine)) {
            trainBests = data.frame(target=factor(breakBestTies(data, i), levels=data$performance))
            if(hasLearnerProperties(combine, "weights") && use.weights) {
                trw = abs(apply(trp, 1, max) - apply(trp, 1, min))
                task = makeClassifTask(id="classifyPairs", target="target", weights=trw, data=data.frame(trainBests, trf$features, trainpredictions), fixup.data="quiet", check.data=FALSE)
            } else {
                task = makeClassifTask(id="classifyPairs", target="target", data=data.frame(trainBests, trf$features, trainpredictions), fixup.data="quiet", check.data=FALSE)
            }
            if(length(unique(trainBests$target)) == 1) {
                # one-class problem
                combinedmodel = train(constantClassifier, task = task)
            } else {
                combinedmodel = train(combine, task = task)
            }
            if(!is.na(save.models)) {
                saveRDS(list(model=combinedmodel, train.data=task, test.data=data.frame(tsf$features, pairpredictions)), file = paste(save.models, combine$id, "combined", i, "rds", sep="."))
            }
            preds = predict(combinedmodel, newdata=data.frame(tsf$features, pairpredictions))
            if(preds$predict.type == "prob") {
                preds = getPredictionProbabilities(preds, preds$task.desc$class.levels)
            } else {
                preds = getPredictionResponse(preds)
                preds = rbind.fill(lapply(preds, function(x) data.frame(t(setNames(as.numeric(x == levels(preds)), levels(preds))))))
            }
            combinedpredictions = rbind.fill(lapply(1:nrow(preds), function(j) {
                ss = preds[j,,drop=F]
                ord = order(ss, decreasing = TRUE)
                data.frame(ids[j,,drop=F], algorithm=names(ss)[ord], score=as.numeric(ss)[ord], iteration=i, row.names = NULL)
            }))
        } else {
            merged = Reduce('+', pairpredictions)
            combinedpredictions = rbind.fill(lapply(1:nrow(merged), function(j) {
                ord = order(merged[j,], decreasing = TRUE)
                data.frame(ids[j,,drop=F], algorithm=names(merged)[ord], score=as.numeric(merged[j,])[ord], iteration=i, row.names = NULL)
            }))
        }
        return(combinedpredictions)
    }, 1:length(data$train), level = "llama.fold"))

    fs = pre(data$data[data$features])
    fp = data$data[data$performance]
    fw = abs(apply(fp, 1, max) - apply(fp, 1, min))
    models = lapply(1:ncol(combns), function(i) {
        if(data$minimize) {
            cmp = function(x, y) {
                sapply(data$data[[x]] < data$data[[y]], function(z) { if(z) { x } else { y } })
            }
        } else {
            cmp = function(x, y) {
                sapply(data$data[[x]] > data$data[[y]], function(z) { if(z) { x } else { y } })
            }
        }
        labels = data.frame(target=factor(cmp(combns[1,i], combns[2,i])))
        if(hasLearnerProperties(classifier, "weights") && use.weights) {
            task = makeClassifTask(id="classifyPairs", target="target", weights=abs(data$data[[combns[1,i]]] - data$data[[combns[2,i]]]), data=data.frame(labels, fs$features), fixup.data="quiet", check.data=FALSE)
        } else {
            task = makeClassifTask(id="classifyPairs", target="target", data=data.frame(labels, fs$features), fixup.data="quiet", check.data=FALSE)
        }
        if(length(unique(labels$target)) == 1) {
            # one-class problem
            model = train(constantClassifier, task = task)
        } else {
            model = train(classifier, task = task)
        }
        return(model)
    })
    if(!is.null(combine)) {
        trainpredictions = list()
        for(i in 1:ncol(combns)) {
            preds = predict(models[[i]], newdata=fs$features)
            trainpredictions[[i]] = if(preds$predict.type == "prob") {
                getPredictionProbabilities(preds, preds$task.desc$class.levels)
            } else {
                tmp = getPredictionResponse(preds)
                rbind.fill(lapply(tmp, function(x) data.frame(t(setNames(as.numeric(x == levels(tmp)), levels(tmp))))))
            }
        }
        if(hasLearnerProperties(combine, "weights") && use.weights) {
            task = makeClassifTask(id="classifyPairs", target="target", weights=fw, data=data.frame(totalBests, fs$features, trainpredictions), fixup.data="quiet", check.data=FALSE)
        } else {
            task = makeClassifTask(id="classifyPairs", target="target", data=data.frame(totalBests, fs$features, trainpredictions), fixup.data="quiet", check.data=FALSE)
        }
        if(length(unique(totalBests$target)) == 1) {
            # one-class problem
            combinedmodel = train(constantClassifier, task = task)
        } else {
            combinedmodel = train(combine, task = task)
        }
    }

    predictor = function(x) {
        tsf = pre(x[data$features], fs$meta)
        if(length(intersect(colnames(x), data$ids)) > 0) {
            ids = x[data$ids]
        } else {
            ids = data.frame(id = 1:nrow(x)) # don't have IDs, generate them
        }
        pairpredictions = list()
        for(i in 1:ncol(combns)) {
            preds = predict(models[[i]], newdata=tsf$features)
            pairpredictions[[i]] = if(preds$predict.type == "prob") {
                getPredictionProbabilities(preds, preds$task.desc$class.levels)
            } else {
                tmp = getPredictionResponse(preds)
                rbind.fill(lapply(tmp, function(x) data.frame(t(setNames(as.numeric(x == levels(tmp)), levels(tmp))))))
            }
        }
        if(!is.null(combine)) {
            preds = predict(combinedmodel, newdata=data.frame(tsf$features, pairpredictions))
            if(preds$predict.type == "prob") {
                preds = getPredictionProbabilities(preds, preds$task.desc$class.levels)
            } else {
                preds = getPredictionResponse(preds)
                preds = rbind.fill(lapply(preds, function(x) data.frame(t(setNames(as.numeric(x == levels(preds)), levels(preds))))))
            }
            combinedpredictions =  rbind.fill(lapply(1:nrow(preds), function(j) {
                ss = preds[j,,drop=F]
                ord = order(ss, decreasing = TRUE)
                data.frame(ids[j,,drop=F], algorithm=names(ss)[ord], score=as.numeric(ss)[ord], iteration=i, row.names = NULL)
            }))
        } else {
            merged = Reduce('+', pairpredictions)
            combinedpredictions = rbind.fill(lapply(1:nrow(merged), function(j) {
                ord = order(merged[j,], decreasing = TRUE)
                data.frame(ids[j,,drop=F], algorithm=names(merged)[ord], score=as.numeric(merged[j,])[ord], iteration=i, row.names = NULL)
            }))
        }
        return(combinedpredictions)
    }
    class(predictor) = "llama.model"
    attr(predictor, "type") = "classifyPairs"
    attr(predictor, "hasPredictions") = FALSE
    attr(predictor, "addCosts") = TRUE

    retval = list(predictions=predictions, models=models, predictor=predictor)
    class(retval) = "llama.model"
    attr(retval, "type") = "classifyPairs"
    attr(retval, "hasPredictions") = TRUE
    attr(retval, "addCosts") = TRUE

    return(retval)
}
class(classifyPairs) = "llama.modelFunction"
