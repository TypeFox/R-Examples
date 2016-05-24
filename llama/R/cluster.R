cluster <-
function(clusterer=NULL, data=NULL, bestBy="performance", pre=function(x, y=NULL) { list(features=x) }, save.models=NA) {
    if(!testClass(clusterer, "Learner") && !testList(clusterer, types="Learner")) {
        stop("Need clusterer or list of clusterers!")
    }
    assertClass(data, "llama.data")
    hs = attr(data, "hasSplits")
    if(is.null(hs) || hs != TRUE) {
        stop("Need data with train/test split!")
    }
    if(testClass(clusterer, "Learner")) { clusterer = list(clusterer) }
    combinator = "majority"
    if(!is.null(clusterer$.combine)) {
        combinator = clusterer$.combine
        clusterer = clusterer[-which(names(clusterer) == ".combine")]
    }

    if(bestBy == "performance") {
        innerbest = function(ss) { sort(sapply(data$performance, function(x) { mean(ss[,x]) }), decreasing=!data$minimize) }
    } else if(bestBy == "count") {
        innerbest = function(ss) { sort(table(data$best[as.integer(rownames(ss))]), decreasing=T) }
    } else if(bestBy == "successes") {
        if(is.null(data$success)) {
            stop("Need successes to determine best by successes!")
        }
        innerbest = function(ss) { setNames(sort(sapply(data$success, function(x) { colMeans(ss[x])[1] }), decreasing=T), data$performance) }
    } else {
        stop(paste("Unknown bestBy: ", bestBy, sep=""))
    }

    worstScore = if(data$minimize) { Inf } else { -Inf }
    bestfun = function(ss) { setNames(data.frame(as.table(innerbest(ss))), predNames) }

    i = 1 # prevent warning when checking package
    predictions = rbind.fill(parallelMap(function(i) {
        trf = pre(data$data[data$train[[i]],][data$features])
        tsf = pre(data$data[data$test[[i]],][data$features], trf$meta)
        ids = data$data[data$test[[i]],][data$ids]

        trainpredictions = matrix(nrow=nrow(trf$features), ncol=length(clusterer))
        ensemblepredictions = list()
        for(j in 1:length(clusterer)) {
            task = makeClusterTask(id="clustering", data=trf$features)
            model = train(clusterer[[j]], task=task)
            if(!is.na(save.models)) {
                saveRDS(list(model=model, train.data=task, test.data=tsf$features), file = paste(save.models, clusterer[[j]]$id, i, "rds", sep="."))
            }
            trainclusters = getPredictionResponse(predict(model, newdata=trf$features))
            if(all(is.na(trainclusters))) {
                best = bestfun(data$data[data$train[[i]],])
            } else {
                best = by(data$data[data$train[[i]],], trainclusters, bestfun)
            }
            if(inherits(combinator, "Learner")) { # only do this if we need it
                trainpredictions[,j] = factor(sapply(getPredictionResponse(predict(model, newdata=trf$features)),
                function(x) {
                    if(is.na(x)) {
                        factor(NA)
                    } else {
                        best[[which(names(best)==x)]][1,1]
                    }
                }))
            }
            preds = getPredictionResponse(predict(model, newdata=tsf$features))
            ensemblepredictions[[j]] = list()
            for(k in 1:length(preds)) {
                x = preds[k]
                if(is.na(x)) {
                    ensemblepredictions[[j]][[k]] = data.frame(algorithm = NA, score = worstScore)
                } else {
                    ensemblepredictions[[j]][[k]] = best[[which(names(best)==x)]]
                }
            }
        }
        if(inherits(combinator, "Learner")) {
            trainBests = data.frame(target=factor(breakBestTies(data, i), levels=data$performance))
            task = makeClassifTask(id="cluster", target="target", data=cbind(trainBests, trf$features, data.frame(trainpredictions)))
            combinedmodel = train(combinator, task=task)
            if(!is.na(save.models)) {
                saveRDS(list(model=combinedmodel, train.data=task, test.data=cbind(tsf$features, data.frame(featureData))), file = paste(save.models, combinator$id, "combined", i, "rds", sep="."))
            }
            featureData = matrix(nrow=nrow(tsf$features), ncol=length(clusterer))
            for(k in 1:nrow(featureData)) {
                for(j in 1:ncol(featureData)) {
                    featureData[k,j] = ensemblepredictions[[j]][[k]][1,1]
                }
            }
            preds = getPredictionResponse(predict(combinedmodel, newdata=cbind(tsf$features, data.frame(featureData))))
            combinedpredictions = rbind.fill(lapply(1:length(preds), function(j) {
                if(all(is.na(preds[j,drop=F]))) {
                    data.frame(ids[j,,drop=F], algorithm=NA, score=worstScore, iteration=i, row.names = NULL)
                } else {
                    tab = table(preds[j,drop=F])
                    data.frame(ids[j,,drop=F], algorithm=names(tab), score=as.vector(tab), iteration=i, row.names = NULL)
                }
            }))
        } else {
            combinedpredictions = rbind.fill(lapply(1:nrow(ids), function(j) {
                df = rbind.fill(lapply(ensemblepredictions, function(x) { x[[j]] }))
                if(all(is.na(df$algorithm))) {
                    data.frame(ids[j,,drop=F], algorithm=NA, score=worstScore, iteration=i, row.names = NULL)
                } else {
                    df = aggregate(data=df, as.formula(paste(predNames[2], predNames[1], sep="~")), sum)
                    df$iteration = i
                    cbind(ids[j,,drop=F], df, row.names = NULL)
                }
            }))
        }
        return(combinedpredictions)
    }, 1:length(data$train), level = "llama.fold"))

    fs = pre(data$data[data$features])
    models = lapply(1:length(clusterer), function(i) {
        task = makeClusterTask(id="clustering", data=fs$features)
        model = train(clusterer[[i]], task=task)
        clusters = getPredictionResponse(predict(model, newdata=fs$features))
        best = by(data$data, clusters, bestfun)
        return(function(newdata) {
            preds = getPredictionResponse(predict(model, newdata=newdata))
            lapply(1:length(preds), function(k) {
                x = preds[k]
                if(is.na(x)) {
                    data.frame(algorithm = NA, score = worstScore)
                } else {
                    best[[which(names(best)==x)]]
                }
            })
        })
    })
    if(inherits(combinator, "Learner")) {
        trainpredictions = matrix(nrow=nrow(fs$features), ncol=length(clusterer))
        for(j in 1:length(clusterer)) {
            trainpredictions[,j] = factor(sapply(models[[j]](fs$features), function(x) { x[1,1] }))
        }
        totalBests = data.frame(target=factor(breakBestTies(data), levels=data$performance))
        combinedmodel = train(combinator, task=makeClassifTask(id="cluster", target="target", data=cbind(totalBests, fs$features, data.frame(trainpredictions))))
    }

    predictor = function(x) {
        tsf = pre(x[data$features], fs$meta)
        if(length(intersect(colnames(x), data$ids)) > 0) {
            ids = x[data$ids]
        } else {
            ids = data.frame(id = 1:nrow(x)) # don't have IDs, generate them
        }
        ensemblepredictions = list()
        for(j in 1:length(clusterer)) {
            ensemblepredictions[[j]] = models[[j]](tsf$features)
        }
        if(inherits(combinator, "Learner")) {
            featureData = matrix(nrow=nrow(tsf$features), ncol=length(clusterer))
            for(k in 1:nrow(featureData)) {
                for(j in 1:ncol(featureData)) {
                    featureData[k,j] = ensemblepredictions[[j]][[k]][1,1]
                }
            }
            preds = getPredictionResponse(predict(combinedmodel, newdata=cbind(tsf$features, data.frame(featureData))))
            combinedpredictions = rbind.fill(lapply(1:length(preds), function(j) {
                if(all(is.na(preds[j,drop=F]))) {
                    data.frame(ids[j,,drop=F], algorithm=NA, score=worstScore, iteration=1, row.names = NULL)
                } else {
                    tab = table(preds[j,drop=F])
                    data.frame(ids[j,,drop=F], algorithm=names(tab), score=as.vector(tab), iteration=1, row.names = NULL)
                }
            }))
        } else {
            combinedpredictions = rbind.fill(lapply(1:nrow(ids), function(j) {
                df = rbind.fill(lapply(ensemblepredictions, function(x) { x[[j]] }))
                if(all(is.na(df$algorithm))) {
                    data.frame(ids[j,,drop=F], algorithm=NA, score=worstScore, iteration=1, row.names = NULL)
                } else {
                    df = aggregate(data=df, as.formula(paste(predNames[2], predNames[1], sep="~")), sum)
                    df$iteration = 1
                    cbind(ids[j,,drop=F], df, row.names = NULL)
                }
            }))
        }
        return(combinedpredictions)
    }
    class(predictor) = "llama.model"
    attr(predictor, "type") = "cluster"
    attr(predictor, "hasPredictions") = FALSE
    attr(predictor, "addCosts") = TRUE

    retval = list(predictions=predictions, models=models, predictor=predictor)
    class(retval) = "llama.model"
    attr(retval, "type") = "cluster"
    attr(retval, "hasPredictions") = TRUE
    attr(retval, "addCosts") = TRUE

    return(retval)
}
class(cluster) = "llama.modelFunction"
