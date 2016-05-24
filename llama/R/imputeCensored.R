imputeCensored <-
function(data=NULL, estimator=makeLearner("regr.lm"), epsilon=0.1, maxit=1000) {
    if(!testClass(estimator, "Learner")) {
        stop("Need regressor to impute values!")
    }
    assertClass(data, "llama.data")
    if(is.null(data$success)) {
        stop("Need successes to impute censored values!")
    }
    if(epsilon <= 0) {
        stop("Epsilon must be > 0!")
    }

    data$original_data = data$data

    i = 0
    for(i in 1:length(data$success)) {
        s = data$success[i]
        p = data$performance[i]
        if(!any(data$data[[s]])) {
            stop(paste("Cannot impute for ", p, ", no non-censored values!"), sep="")
        }
        if(!all(data$data[[s]])) {
            splits = split(1:nrow(data$data), data$data[[s]])
            haveind = splits$`TRUE`
            wantind = splits$`FALSE`
            task = makeRegrTask(id="imputation", target="target", data=cbind(data.frame(target=data$data[haveind,p]), data$data[haveind,][data$features]))
            model = train(estimator, task=task)
            data$data[wantind,p] = predict(model, newdata=data$data[wantind,][data$features])$data$response

            diff = Inf
            it = 1
            while(diff > epsilon) {
                task = makeRegrTask(id="imputation", target="target", data=cbind(data.frame(target=data$data[[p]]), data$data[data$features]))
                model = train(estimator, task=task)
                preds = predict(model, newdata=data$data[wantind,][data$features])$data$response
                diff = max(abs(preds - data$data[wantind,p]))
                data$data[wantind,p] = preds
                it = it + 1
                if(it > maxit) {
                    warning(paste("Did not reach convergence within ", maxit, " iterations for ", p, ".", sep=""))
                    break
                }
            }
            data$data[[s]] = rep.int(T, nrow(data$data))
        }
    }

    return(data)
}
