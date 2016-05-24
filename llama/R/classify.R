classify <-
function(classifier=NULL, data=NULL, pre=function(x, y=NULL) { list(features=x) }, save.models=NA, use.weights = TRUE) {
    if(!testClass(classifier, "Learner") && !testList(classifier, types="Learner")) {
        stop("Need classifier or list of classifiers!")
    }
    assertClass(data, "llama.data")
    hs = attr(data, "hasSplits")
    if(is.null(hs) || hs != TRUE) {
        stop("Need data with train/test split!")
    }
    if(testClass(classifier, "Learner")) { classifier = list(classifier) }
    combinator = "majority"
    if(!is.null(classifier$.combine)) {
        combinator = classifier$.combine
        classifier = classifier[-which(names(classifier) == ".combine")]
    }

    totalBests = data.frame(target=factor(breakBestTies(data), levels=data$performance))

    predictions = rbind.fill(parallelMap(function(i) {
        trf = pre(data$data[data$train[[i]],][data$features])
        tsf = pre(data$data[data$test[[i]],][data$features], trf$meta)
        ids = data$data[data$test[[i]],][data$ids]
        trp = data$data[data$train[[i]],][data$performance]
        trw = abs(apply(trp, 1, max) - apply(trp, 1, min))

        trainpredictions = list()
        ensemblepredictions = list()

        trainBests = data.frame(target=factor(breakBestTies(data, i), levels=data$performance))
        for(j in 1:length(classifier)) {
            if(hasLearnerProperties(classifier[[j]], "weights") && use.weights) {
                task = makeClassifTask(id="classify", target="target", weights = trw, data=data.frame(trainBests, trf$features))
            } else {
                task = makeClassifTask(id="classify", target="target", data=data.frame(trainBests, trf$features))
            }
            if(length(unique(trainBests$target)) == 1) {
                # one-class problem
                model = train(constantClassifier, task = task)
            } else {
                model = train(classifier[[j]], task = task)
            }
            if(!is.na(save.models)) {
                saveRDS(list(model=model, train.data=task, test.data=tsf$features), file = paste(save.models, classifier[[j]]$id, i, "rds", sep="."))
            }
            if(inherits(combinator, "Learner")) { # only do this if we need it
                preds = predict(model, newdata=trf$features)
                trainpredictions[[j]] = if(preds$predict.type == "prob") {
                    getPredictionProbabilities(preds, data$performance)
                } else {
                    tmp = getPredictionResponse(preds)
                    rbind.fill(lapply(tmp, function(x) data.frame(t(setNames(as.numeric(x == levels(tmp)), levels(tmp))))))
                }
            }
            preds = predict(model, newdata=tsf$features)
            ensemblepredictions[[j]] = if(preds$predict.type == "prob") {
                getPredictionProbabilities(preds, data$performance)
            } else {
                tmp = getPredictionResponse(preds)
                rbind.fill(lapply(tmp, function(x) data.frame(t(setNames(as.numeric(x == levels(tmp)), levels(tmp))))))
            }
        }
        if(inherits(combinator, "Learner")) {
            if(hasLearnerProperties(combinator, "weights") && use.weights) {
                task = makeClassifTask(id="classify", target="target", weights = trw, data=data.frame(trainBests, trf$features, trainpredictions))
            } else {
                task = makeClassifTask(id="classify", target="target", data=data.frame(trainBests, trf$features, trainpredictions))
            }
            if(length(unique(trainBests$target)) == 1) {
                # one-class problem
                combinedmodel = train(constantClassifier, task = task)
            } else {
                combinedmodel = train(combinator, task = task)
            }
            if(!is.na(save.models)) {
                saveRDS(list(model=combinedmodel, train.data=task, test.data=data.frame(tsf$features, ensemblepredictions)), file = paste(save.models, combinator$id, "combined", i, "rds", sep="."))
            }
            preds = predict(combinedmodel, newdata=data.frame(tsf$features, ensemblepredictions))
            if(preds$predict.type == "prob") {
                preds = getPredictionProbabilities(preds, data$performance)
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
            merged = Reduce('+', ensemblepredictions)
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
    models = lapply(1:length(classifier), function(i) {
        if(hasLearnerProperties(classifier[[i]], "weights") && use.weights) {
            task = makeClassifTask(id="classify", target="target", weights = fw, data=data.frame(totalBests, fs$features))
        } else {
            task = makeClassifTask(id="classify", target="target", data=data.frame(totalBests, fs$features))
        }
        if(length(unique(totalBests$target)) == 1) {
            # one-class problem
            model = train(constantClassifier, task = task)
        } else {
            model = train(classifier[[i]], task = task)
        }
        return(model)
    })
    if(inherits(combinator, "Learner")) {
        trainpredictions = list()
        for(i in 1:length(classifier)) {
            preds = predict(models[[i]], newdata=fs$features)
            trainpredictions[[i]] = if(preds$predict.type == "prob") {
                getPredictionProbabilities(preds, data$performance)
            } else {
                tmp = getPredictionResponse(preds)
                rbind.fill(lapply(tmp, function(x) data.frame(t(setNames(as.numeric(x == levels(tmp)), levels(tmp))))))
            }
        }
        if(hasLearnerProperties(combinator, "weights") && use.weights) {
            task = makeClassifTask(id="classify", target="target", weights = fw, data=data.frame(totalBests, fs$features, trainpredictions))
        } else {
            task = makeClassifTask(id="classify", target="target", data=data.frame(totalBests, fs$features, trainpredictions))
        }
        if(length(unique(totalBests$target)) == 1) {
            # one-class problem
            combinedmodel = train(constantClassifier, task = task)
        } else {
            combinedmodel = train(combinator, task = task)
        }
    }

    predictor = function(x) {
        tsf = pre(x[data$features], fs$meta)
        if(length(intersect(colnames(x), data$ids)) > 0) {
            ids = x[data$ids]
        } else {
            ids = data.frame(id = 1:nrow(x)) # don't have IDs, generate them
        }
        ensemblepredictions = list()
        for(i in 1:length(classifier)) {
            preds = predict(models[[i]], newdata=tsf$features)
            ensemblepredictions[[i]] = if(preds$predict.type == "prob") {
                tmp = getPredictionProbabilities(preds, data$performance)
            } else {
                tmp = getPredictionResponse(preds)
                rbind.fill(lapply(tmp, function(x) data.frame(t(setNames(as.numeric(x == levels(tmp)), levels(tmp))))))
            }
        }
        if(inherits(combinator, "Learner")) {
            preds = predict(combinedmodel, newdata=data.frame(tsf$features, ensemblepredictions))
            if(preds$predict.type == "prob") {
                preds = getPredictionProbabilities(preds, data$performance)
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
            merged = Reduce('+', ensemblepredictions)
            combinedpredictions = rbind.fill(lapply(1:nrow(merged), function(j) {
                ord = order(merged[j,], decreasing = TRUE)
                data.frame(ids[j,,drop=F], algorithm=names(merged)[ord], score=as.numeric(merged[j,])[ord], iteration=i, row.names = NULL)
            }))
        }
        return(combinedpredictions)
    }
    class(predictor) = "llama.model"
    attr(predictor, "type") = "classify"
    attr(predictor, "hasPredictions") = FALSE
    attr(predictor, "addCosts") = TRUE

    retval = list(predictions=predictions, models=models, predictor=predictor)
    class(retval) = "llama.model"
    attr(retval, "type") = "classify"
    attr(retval, "hasPredictions") = TRUE
    attr(retval, "addCosts") = TRUE

    return(retval)
}
class(classify) = "llama.modelFunction"
